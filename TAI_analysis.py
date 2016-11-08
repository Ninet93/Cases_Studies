import pandas as pd
import numpy as np
import GEOparse
import matplotlib.pyplot as plt

probes_conv = GEOparse.parse_GSM("/Users/Baboo19/Documents/Master/Automne_16/Case_Studies/module1/data/GPL6457_old_annotations.txt.gz")

gse = GEOparse.get_GEO("GSE24616", destdir="./")


char = {"stage": [], "time": [], "sex": [], "sample_name": []}
for gsm_name, gsm in sorted(gse.gsms.iteritems()):
    char["stage"].append(gsm.metadata['characteristics_ch1'][1].split(": ")[1])
    char["time"].append(gsm.metadata['characteristics_ch1'][2].split(": ")[1])
    char["sex"].append(gsm.metadata['characteristics_ch1'][3].split(": ")[1])
    char["sample_name"].append(gsm.name)


print(char["stage"][0:10], char["time"][0:10],char["sex"][0:10], char["sample_name"][0:10])

gsm.head()


GPL = gse.gpls["GPL6457"]
pivoted_samples = gse.pivot_samples('VALUE')
#find matching genes in GPL transcripts



pivoted_samples.set_index(GPL.table.SPOT_ID, inplace=True)

# read in age index
strata = pd.read_csv("/Users/Baboo19/Documents/Master/Automne_16/Case_Studies/module1/data/processed/danio_age_index.txt",sep="\t",header=None)
strata.columns = ["GeneID","ProbeID","age"]
strata.set_index("ProbeID",inplace=True)
# age = ps

#find genes in AI file that match expression data ('inner', refers to intersect)
matched_data = pivoted_samples.join(strata, how="inner").groupby(level=0).last()


#average out multiple transcripts
unique_data = matched_data.groupby("GeneID").mean()


#only used mixed and female samples
char_pd = pd.DataFrame(char,index=char["sample_name"])
mixed = char_pd[char_pd.sex == "mixed"].sample_name.tolist()
mixed += char_pd[char_pd.sex == "female"].sample_name.tolist()


#sort according time points
char_pd["timing_number"] = 0
time_stamps = char_pd.time.unique()
for i in xrange(len(time_stamps)):
    char_pd.loc[char_pd.time == time_stamps[i],"timing_number"] = i+1


#average the samples across identicla time points
experiment_index = char_pd[char_pd.index.isin(mixed)].reset_index().groupby("timing_number")["index"].apply(lambda x: np.array(x))
set_mean = {}
stages = []
for d, col_list in experiment_index.iteritems():
    set_mean[d] = unique_data[col_list].mean(axis=1)
    stages.append(char_pd[char_pd.index.isin(col_list)].stage[0])

mean_data = pd.DataFrame(set_mean)
mean_data['age'] = unique_data['age']
# equivalent to the processed data ProcessedMicroarrayData.txt - almost !

mean_data.head()
mean_data.info()







# TAI = (ps1 * E1 + ps2 * E2) / E1 * E2


unique_data.iloc[1,3] # row, column



TAI=[]
for j in range(0,len(mean_data.columns)-1): # columns # 60
    num=0
    den=0
    for i in range(0,len(mean_data['age'])): # rows # 12891
        num=mean_data['age'][i]*mean_data.iloc[i,j] + num
        den=den + mean_data.iloc[i,j]
        

    TAI.append(num/den)



plt.plot(TAI)
plt.show()
plt.savefig('TAI.png')



del mean_data['age']


mean_data.hist()



mat = mean_data.as_matrix() # not used


val=[]
for j in range(0, len(mean_data.columns)):
    for i in range(0, len(mean_data[1])):
        val.append(mean_data.iloc[i,j])
   
        
plt.hist(val)
plt.hist(val, bins=100)


max(val)

np.log(val)





#### LOG TRANSFORMATION OF THE DATA

del mean_data['age']

log_TAI=[]
for j in range(0,len(mean_data.columns)-1): # columns # 60
    num=0
    den=0
    for i in range(0,len(mean_data['age'])): # rows # 12891
        num=mean_data['age'][i]*np.log(mean_data.iloc[i,j]) + num
        den=den + np.log(mean_data.iloc[i,j])
        

    log_TAI.append(num/den)
# to tranform the data with log() isn't the same as doing an histogram and put log=True !!!!!

plt.plot(log_TAI)





log_val=[]
for j in range(0, len(mean_data.columns)):
    for i in range(0, len(mean_data[1])):
        log_val.append(np.log(mean_data.iloc[i,j]))
   
        
plt.hist(log_val)
plt.hist(log_val, bins=100)

