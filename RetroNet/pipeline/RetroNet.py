import os
import sys
import shutil
import numpy as np
import torch
import torch.optim as optim
import torchvision
import time, datetime
import copy
import random
from torchvision.models import resnet18
from torch.utils.data import Dataset, DataLoader
from torch import nn
from PIL import Image

# Get parameters #
outpath=sys.argv[1]
sub=sys.argv[2]
TEclass=sys.argv[3]
ver=sys.argv[4]
masterpath=sys.argv[5]
cutoff=sys.argv[6]
hg=sys.argv[7]

### PNG list ###
Samples_dir = outpath + '/' + sub + '/visual_' + ver + '/'+ TEclass
sample_list = []
for dirpath, dirnames, filenames in os.walk(Samples_dir):
    for filename in filenames:
        if Image.open(os.path.join(dirpath, filename)).convert('RGB').size == (6620,60) and os.path.getsize(os.path.join(dirpath, filename)) != 0:  
            sample_list.append(os.path.join(dirpath, filename))
print(len(sample_list),TEclass,"images founded")
### Custom my own dataset ###
class MyDataset(Dataset):
    def __init__(self, file_list):
        samples_list = []
        samples_label = []
        for sample_file in file_list:
            samples_list.append(sample_file)
            if "TRUE" in sample_file:
                samples_label.append(1)
            else:
                samples_label.append(0)
        self.data = np.array(samples_list)
        self.labels = np.array(samples_label)
        self.tensor = torchvision.transforms.Compose([torchvision.transforms.ToTensor()])
    def __getitem__(self, index):
        img = self.data[index]
        label = self.labels[index]
        img = Image.open(img).convert('RGB')
        data = self.tensor(img)
        return data, label
    def __len__(self): #necessary function
        return self.data.shape[0]

## Choose model ##
if TEclass == "LINE":
    Best_model_dir = masterpath + "/RetroNet/LINE_RetroNet_model.pt"
    Inspect_script = "./LINE_Inpection.py"
elif TEclass == "ALU":
    Best_model_dir = masterpath + "/RetroNet/ALU_RetroNet_model.pt"
    Inspect_script = "./ALU_Inpection.py"
elif TEclass == "SVA":
    Best_model_dir = masterpath + "/RetroNet/SVA_RetroNet_model.pt"
    Inspect_script = "./SVA_Inpection.py"
### Initial Setting ###
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
net = resnet18(weights=None, progress=True,num_classes=2)
net.load_state_dict(torch.load(Best_model_dir,map_location=torch.device(device)))
net.to(device)
net.eval()

### Prediction ###
test_dataset = MyDataset(sample_list)
test_dataloader = DataLoader(test_dataset,batch_size=24,shuffle=False)
y_pred = []
y_scores = []
for i, (inputs, labels) in enumerate(test_dataloader):
        inputs = inputs.to(device)
        preds = net(inputs)
        probs = nn.functional.softmax(preds,dim=1)
        _,predictions = torch.max(probs,1)
        y_scores.extend(probs.detach().cpu().numpy()[:,1])
        y_pred.extend(predictions.cpu().numpy())
# output all probability #
f = open("{}.txt".format(outpath + "/" + sub + "/RetroNet/" + TEclass + "_Probability"),"w")
f.write("Image\tProbability\n")
for i in range(len(y_scores)):
    f.write("{}\t{:.10f}\n".format(sample_list[i],y_scores[i]))
f.close()

# Inspect probs larger than cutoff #
f = open("{}.txt".format(outpath + "/" + sub + "/RetroNet/" + TEclass + "_Probability"),"r")
records = f.readlines()
f = open("{}_insertion_site.bed".format(TEclass),"w")
for line in records[1:]:
    f.write("{}\t{}\t{}\t{}\n".format(line.split("\n")[0].split("/")[-1].split("_")[-3],line.split("\n")[0].split("/")[-1].split("_")[-2],int(line.split("\n")[0].split("/")[-1].split("_")[-2])+1,line.split("\n")[0]))
f.close()
if hg == "hg38":
    os.system("windowBed -a {}_insertion_site.bed -b {}/RetroNet/sorted.hg38.centr.telos.pos.bed -w 1000 -c > {}_insertion_site_hit.bed".format(TEclass,masterpath,TEclass))
elif hg == "b37":
    os.system("windowBed -a {}_insertion_site.bed -b {}/RetroNet/sorted.b37.centr.telos.pos.bed -w 1000 -c > {}_insertion_site_hit.bed".format(TEclass,masterpath,TEclass))
os.system("awk '{}' {}_insertion_site_hit.bed > {}_insertion_site_hited.bed".format("{print $4, $5, $6, $7}",TEclass,TEclass))
os.system("rm -f {}_insertion_site.bed {}_insertion_site_hit.bed".format(TEclass,TEclass))

f = open("{}_insertion_site_hited.bed".format(TEclass),"r")
filtered_record = f.readlines()
os.system("rm {}_insertion_site_hited.bed".format(TEclass))

f = open("{}/{}_Inspected_{}_cut{}.txt".format(outpath + "/" + sub + "/RetroNet",TEclass,sub,cutoff),"w")
for line in filtered_record:
    if float(line.split()[1]) >= float(cutoff):
        marker = ""
        try:
            os.system("python3 {} {} {}".format(Inspect_script,line.split()[0],ver))
            f1 = open("Inspection_tmp_{}_{}.txt".format(TEclass,ver),"r")
            inspection = f1.readline()
        except:
            inspection = "Can't perform inspection!!!"
        line_split = line.split()
        if int(line_split[2]) != 0:
            marker = "Centromere or Telemere"
        else:
            marker = "0"
        f.write(line_split[0] + "\t" + line_split[1] + "\t" + marker + "\t" + inspection + "\n")
f.close()
