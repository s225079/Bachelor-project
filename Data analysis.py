


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.scale import ScaleBase
from matplotlib.transforms import Transform
from matplotlib import scale as mscale
from matplotlib.ticker import MultipleLocator

import matplotlib.cm as cm
import matplotlib.colors as mcolors

#Substract data from excel files

file_path_t0 = 'C:\Od data\S.cop_210525_t0.xlsx'
file_path_t24 = 'C:\Od data\s.cop_220525_t24.xlsx'
file_path_t48 = 'C:\Od data\S.cop_230525_t48.xlsx'
# Read specific Excel range: C29:L34
df_t0 = pd.read_excel(
    file_path_t0,
    sheet_name='Sheet2',   # Change if needed
    engine='openpyxl',
    usecols='C:L',
    skiprows=27,    # Skip rows before row 29
    nrows=6         # Rows 29–34
)

df_t24 = pd.read_excel(
    file_path_t24,
    sheet_name='Sheet2',   # Change if needed
    engine='openpyxl',
    usecols='C:L',
    skiprows=27,    # Skip rows before row 29
    nrows=6         # Rows 29–34
)
print(df_t24)
df_t48 = pd.read_excel(
    file_path_t48,
    sheet_name='Sheet2',   # Change if needed
    engine='openpyxl',
    usecols='C:L',
    skiprows=27,    # Skip rows before row 29
    nrows=6         # Rows 29–34
)

# Convert to NumPy array
array_t0 = df_t0.to_numpy()
array_t24 = df_t24.to_numpy()
array_t48 = df_t48.to_numpy()

#Put data into a vector:

data_t0=np.array([])
data_t24=np.array([])
data_t48=np.array([])


for i in range (0,6):
    
    data_t0=np.concatenate((data_t0,(array_t0[i,1:5]-array_t0[i,0])))
    data_t24=np.concatenate((data_t24,(array_t24[i,1:5]-array_t24[i,0])))
    data_t48=np.concatenate((data_t48,(array_t48[i,1:5]-array_t48[i,0])))

    
for i in range (0,6):
    
    data_t0=np.concatenate((data_t0,(array_t0[i,5:9]-array_t0[i,9])))
    data_t24=np.concatenate((data_t24,(array_t24[i,5:9]-array_t24[i,9])))
    data_t48=np.concatenate((data_t48,(array_t48[i,5:9]-array_t48[i,9])))
    


#Merge vectors to create the whole dataset with all timepoints:

all_data=np.vstack((data_t0,data_t24,data_t48))
all_data=np.transpose(all_data)

  
#Make plots with the mean and std for each timeframe:

mean_H20=np.array([])
std_H20=np.array([])

mean_Glu=np.array([])
std_Glu=np.array([])

mean_Ara=np.array([])
std_Ara=np.array([])

mean_SBP=np.array([])
std_SBP=np.array([])

mean_CP=np.array([])
std_CP=np.array([])

mean_RGI=np.array([])
std_RGI=np.array([])

mean_GalMan=np.array([])
std_GalMan=np.array([])

mean_FOS=np.array([])
std_FOS=np.array([])

mean_Inulin=np.array([])
std_Inulin=np.array([])

mean_Starch=np.array([])
std_Starch=np.array([])

mean_bglucan=np.array([])
std_bglucan=np.array([])

mean_H20_2=np.array([])
std_H20_2=np.array([])

for i in range (0,3):
    #water
    mean_H20=np.concatenate((mean_H20,np.array([np.mean(all_data[0:4,i])])))
    std_H20=np.concatenate((std_H20,np.array([np.std(all_data[0:4,i])])))

    #Gluc
    mean_Glu=np.concatenate((mean_Glu,np.array([np.mean(all_data[4:8,i])])))
    std_Glu=np.concatenate((std_Glu,np.array([np.std(all_data[4:8,i])])))
    #Ara
    mean_Ara=np.concatenate((mean_Ara,np.array([np.mean(all_data[8:12,i])])))
    std_Ara=np.concatenate((std_Ara,np.array([np.std(all_data[8:12,i])])))
    #SBP
    mean_SBP=np.concatenate((mean_SBP,np.array([np.mean(all_data[12:16,i])])))
    std_SBP=np.concatenate((std_SBP,np.array([np.std(all_data[12:16,i])])))
    #CP
    mean_CP=np.concatenate((mean_CP,np.array([np.mean(all_data[16:20,i])])))
    std_CP=np.concatenate((std_CP,np.array([np.std(all_data[16:20,i])])))
    #RG-I
    mean_RGI=np.concatenate((mean_RGI,np.array([np.mean(all_data[20:24,i])])))
    std_RGI=np.concatenate((std_RGI,np.array([np.std(all_data[20:24,i])])))
    #GalMan
    mean_GalMan=np.concatenate((mean_GalMan,np.array([np.mean(all_data[24:28,i])])))
    std_GalMan=np.concatenate((std_GalMan,np.array([np.std(all_data[24:28,i])])))
    #FOS
    mean_FOS=np.concatenate((mean_FOS,np.array([np.mean(all_data[28:32,i])])))
    std_FOS=np.concatenate((std_FOS,np.array([np.std(all_data[28:32,i])])))
    #Inulin
    mean_Inulin=np.concatenate((mean_Inulin,np.array([np.mean(all_data[32:36,i])])))
    std_Inulin=np.concatenate((std_Inulin,np.array([np.std(all_data[32:36,i])])))
    #Starch
    mean_Starch=np.concatenate((mean_Starch,np.array([np.mean(all_data[36:40,i])])))
    std_Starch=np.concatenate((std_Starch,np.array([np.std(all_data[36:40,i])])))
    #beta-glucan
    mean_bglucan=np.concatenate((mean_bglucan,np.array([np.mean(all_data[40:44,i])])))
    std_bglucan=np.concatenate((std_bglucan,np.array([np.std(all_data[40:44,i])])))
    #water(control)
    mean_H20_2=np.concatenate((mean_H20_2,np.array([np.mean(all_data[44:48,i])])))
    std_H20_2=np.concatenate((std_H20_2,np.array([np.std(all_data[44:48,i])])))

    

t=[0,24,48]

# Plot

plt.errorbar(t, mean_Glu, yerr=std_Glu, fmt='o-', capsize=5, label='Glucose', linewidth=3)
plt.errorbar(t, mean_H20, yerr=std_H20, fmt='o-', capsize=5, label='Water',linewidth=3)
plt.errorbar(t, mean_Ara, yerr=std_Ara, fmt='o-', capsize=5, label='Arabino-oligosaccharides',linewidth=3)
plt.errorbar(t, mean_SBP, yerr=std_SBP, fmt='o-', capsize=5, label='Sugar beet pectin',linewidth=3)
plt.errorbar(t, mean_CP, yerr=std_CP, fmt='o-', capsize=5, label='Citrus pectin',linewidth=3)
plt.errorbar(t, mean_RGI, yerr=std_RGI, fmt='o-', capsize=5, label='RG-I',linewidth=3)
plt.errorbar(t, mean_GalMan, yerr=std_GalMan, fmt='o-', capsize=5, label='Galactomannan',linewidth=3)
plt.errorbar(t, mean_FOS, yerr=std_FOS, fmt='o-', capsize=5, label='Fructo-oligosaccharides',linewidth=3)
plt.errorbar(t, mean_Inulin, yerr=std_Inulin, fmt='o-', capsize=5, label='Inulin',linewidth=3)
plt.errorbar(t, mean_Starch, yerr=std_Starch, fmt='o-', capsize=5, label='Starch',linewidth=3)
plt.errorbar(t, mean_bglucan, yerr=std_bglucan, fmt='o-', capsize=5, label='beta-glucan', color="blue",linewidth=3)

# Labels and title
plt.xlabel('Time (hours)', fontsize=20)
plt.ylabel('Absorbance (OD_600)', fontsize=20)
plt.title('Absorbance measurements over time (S. copri)', fontsize=24)
plt.grid(True)
plt.xticks(t, fontsize=18)
plt.yticks(fontsize=18)

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=20, title="Carbon source", title_fontsize=20)
plt.subplots_adjust(left=0.1, bottom=0.09, right=0.686, top=0.9)

plt.show()



######################## Relatiev change of OD #########################
#Put data into a vector:

data_t0=np.array([])
data_t24=np.array([])
data_t48=np.array([])

for i in range (0,6):
    
    data_t0=np.concatenate((data_t0,(array_t0[i,1:5]-array_t0[i,0])))
    #data_t6=np.concatenate((data_t6,(array_t6[i,1:5]-array_t6[i,0])))
    data_t24=np.concatenate((data_t24,(array_t24[i,1:5]-array_t24[i,0])))
    data_t48=np.concatenate((data_t48,(array_t48[i,1:5]-array_t48[i,0])))

    
for i in range (0,6):
    
    data_t0=np.concatenate((data_t0,(array_t0[i,5:9]-array_t0[i,9])))
    #data_t6=np.concatenate((data_t6,(array_t6[i,5:9]-array_t6[i,9])))
    data_t24=np.concatenate((data_t24,(array_t24[i,5:9]-array_t24[i,9])))
    data_t48=np.concatenate((data_t48,(array_t48[i,5:9]-array_t48[i,9])))
    

#Calculate relative change




data_OD_t24 = (data_t24 - data_t0)/np.absolute(data_t0)
data_OD_t48 = (data_t48 - data_t0)/np.absolute(data_t0)
    


#Merge vectors to create the whole dataset with all timepoints:

all_OD_data=np.vstack((data_OD_t24,data_OD_t48))
all_OD_data=np.transpose(all_OD_data)
np.shape(all_OD_data)
all_OD_data=np.hstack((np.zeros((48,1)),all_OD_data))



#Make plots with the mean and std for each timeframe:

mean_H20=np.array([])
std_H20=np.array([])

mean_Glu=np.array([])
std_Glu=np.array([])

mean_Ara=np.array([])
std_Ara=np.array([])

mean_SBP=np.array([])
std_SBP=np.array([])

mean_CP=np.array([])
std_CP=np.array([])

mean_RGI=np.array([])
std_RGI=np.array([])

mean_GalMan=np.array([])
std_GalMan=np.array([])

mean_FOS=np.array([])
std_FOS=np.array([])

mean_Inulin=np.array([])
std_Inulin=np.array([])

mean_Starch=np.array([])
std_Starch=np.array([])

mean_bglucan=np.array([])
std_bglucan=np.array([])

mean_H20_2=np.array([])
std_H20_2=np.array([])



for i in range (1,3):
    #water
    mean_H20 = np.concatenate((mean_H20, np.array([np.mean(all_OD_data[0:4, i])])))
    std_H20 = np.concatenate((std_H20, np.array([np.std(all_OD_data[0:4, i])])))
    # Gluc
    mean_Glu = np.concatenate((mean_Glu, np.array([np.mean(all_OD_data[4:8, i])])))
    std_Glu = np.concatenate((std_Glu, np.array([np.std(all_OD_data[4:8, i])])))
    # Ara
    mean_Ara = np.concatenate((mean_Ara, np.array([np.mean(all_OD_data[8:12, i])])))
    std_Ara = np.concatenate((std_Ara, np.array([np.std(all_OD_data[8:12, i])])))
    # SBP
    mean_SBP = np.concatenate((mean_SBP, np.array([np.mean(all_OD_data[12:16, i])])))
    std_SBP = np.concatenate((std_SBP, np.array([np.std(all_OD_data[12:16, i])])))
    # CP
    mean_CP = np.concatenate((mean_CP, np.array([np.mean(all_OD_data[16:20, i])])))
    std_CP = np.concatenate((std_CP, np.array([np.std(all_OD_data[16:20, i])])))
    # RG-I
    mean_RGI = np.concatenate((mean_RGI, np.array([np.mean(all_OD_data[20:24, i])])))
    std_RGI = np.concatenate((std_RGI, np.array([np.std(all_OD_data[20:24, i])])))
    # GalMan
    mean_GalMan = np.concatenate((mean_GalMan, np.array([np.mean(all_OD_data[24:28, i])])))
    std_GalMan = np.concatenate((std_GalMan, np.array([np.std(all_OD_data[24:28, i])])))
    # FOS
    mean_FOS = np.concatenate((mean_FOS, np.array([np.mean(all_OD_data[28:32, i])])))
    std_FOS = np.concatenate((std_FOS, np.array([np.std(all_OD_data[28:32, i])])))
    # Inulin
    mean_Inulin = np.concatenate((mean_Inulin, np.array([np.mean(all_OD_data[32:36, i])])))
    std_Inulin = np.concatenate((std_Inulin, np.array([np.std(all_OD_data[32:36, i])])))
    # Starch
    mean_Starch = np.concatenate((mean_Starch, np.array([np.mean(all_OD_data[36:40, i])])))
    std_Starch = np.concatenate((std_Starch, np.array([np.std(all_OD_data[36:40, i])])))
    # beta-glucan
    mean_bglucan = np.concatenate((mean_bglucan, np.array([np.mean(all_OD_data[40:44, i])])))
    std_bglucan = np.concatenate((std_bglucan, np.array([np.std(all_OD_data[40:44, i])])))
    # water (control)
    mean_H20_2 = np.concatenate((mean_H20_2, np.array([np.mean(all_OD_data[44:48, i])])))
    std_H20_2 = np.concatenate((std_H20_2, np.array([np.std(all_OD_data[44:48, i])])))





# --- Define custom scale class ---
class SignedSqrtScale(ScaleBase):
    name = 'signed_sqrt'

    def get_transform(self):
        return self.SignedSqrtTransform()

    def set_default_locators_and_formatters(self, axis):
        pass

    class SignedSqrtTransform(Transform):
        input_dims = output_dims = 1

        def transform_non_affine(self, a):
            return np.sign(a) * np.sqrt(np.abs(a))

        def inverted(self):
            return SignedSqrtScale.InvertedSignedSqrtTransform()

    class InvertedSignedSqrtTransform(Transform):
        input_dims = output_dims = 1

        def transform_non_affine(self, a):
            return np.sign(a) * (a**2)

        def inverted(self):
            return SignedSqrtScale.SignedSqrtTransform()



mscale.register_scale(SignedSqrtScale)



t=[24,48]


# Plot

plt.errorbar(t, mean_Glu, yerr=std_Glu, fmt='o-', capsize=10, elinewidth=3, markersize=5, linewidth=3, alpha=0.9, capthick=3, label='Glucose')
plt.errorbar(t, mean_H20, yerr=std_H20, fmt='o-', capsize=10, elinewidth=3, markersize=5, linewidth=3, alpha=0.9, capthick=3, label='Water')
plt.errorbar(t, mean_Ara, yerr=std_Ara, fmt='o-', capsize=10, elinewidth=3, markersize=5, linewidth=3, alpha=0.9, capthick=3, label='Arabino-oligosaccharides')
plt.errorbar(t, mean_SBP, yerr=std_SBP, fmt='o-', capsize=10, elinewidth=3, markersize=5, linewidth=3, alpha=0.9,capthick=3, label='Sugar beet pectin')
plt.errorbar(t, mean_CP, yerr=std_CP, fmt='o-', capsize=10, elinewidth=3, markersize=5, linewidth=3, alpha=0.9,capthick=3, label='Citrus pectin')
plt.errorbar(t, mean_RGI, yerr=std_RGI, fmt='o-', capsize=10, elinewidth=3, markersize=5, linewidth=3, alpha=0.9, capthick=3,label='RG-I')
plt.errorbar(t, mean_GalMan, yerr=std_GalMan, fmt='o-', capsize=10, elinewidth=3, markersize=5, linewidth=3, alpha=0.9,capthick=3, label='Galactomannan')
plt.errorbar(t, mean_FOS, yerr=std_FOS, fmt='o-', capsize=10, elinewidth=3, markersize=5, linewidth=3, alpha=0.9, capthick=3,label='Fructo-oligosaccharides')
plt.errorbar(t, mean_Inulin, yerr=std_Inulin, fmt='o-', capsize=10, elinewidth=3, markersize=5, linewidth=3, alpha=0.9,capthick=3, label='Inulin')
plt.errorbar(t, mean_Starch, yerr=std_Starch, fmt='o-', capsize=10, elinewidth=3, markersize=5, linewidth=3, alpha=0.9,capthick=3, label='Starch')
plt.errorbar(t, mean_bglucan, yerr=std_bglucan, fmt='o-', capsize=10, elinewidth=3, markersize=5, linewidth=3, alpha=0.9, capthick=3,label='beta-glucan', color="blue")












# Your existing data placeholders and variables
x = np.array([0, 1.5])  # Positions for 24h and 48h
t = ['24 h', '48 h']
bar_width = 0.08

means = [mean_H20, mean_Glu, mean_Ara, mean_SBP, mean_CP, mean_RGI,
         mean_GalMan, mean_FOS, mean_Inulin, mean_Starch, mean_bglucan]

stds = [std_H20, std_Glu, std_Ara, std_SBP, std_CP, std_RGI,
        std_GalMan, std_FOS, std_Inulin, std_Starch, std_bglucan]

labels = [ 'Water', 'Glucose', 'AOS', 'SBP', 'CP', 'RG-I',
          'Galactomannan', 'FOS', 'Inulin', 'Starch', 'beta-glucan']

# Define colors (you can choose any matplotlib color names or hex codes)
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
          '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', "#1831BE"]

# Plot bars with assigned colors
for i in range(len(means)):
    offset = (i - 5) * bar_width
    plt.bar(x + offset, means[i], width=bar_width, yerr=stds[i],
            capsize=7, label=labels[i], error_kw=dict(lw=3.5), color=colors[i])

# Formatting (same as before)
plt.xticks(x, t)
plt.xlabel('Time (hours)', fontsize=26)
plt.ylabel('Relative ΔOD_600', fontsize=26)
plt.title('Relative change in absorbance measurements over time ($\\mathit{S.\\ copri}$)', fontsize=28)
plt.grid(True, axis='y')
plt.xticks(fontsize=24)
plt.yticks(fontsize=22)
plt.axhline(0, color='black', linewidth=2)
plt.yscale('signed_sqrt')
plt.gca().yaxis.set_major_locator(MultipleLocator(1))

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=24, title="Carbon source", title_fontsize=24)
plt.grid(visible=True, which='major', axis='y', linestyle='--', linewidth=1.5)
plt.grid(visible=True, which='minor', axis='y', linestyle='--', linewidth=1.5)
plt.subplots_adjust(left=0.075, bottom=0.1, right=0.78, top=0.9)
plt.show()









