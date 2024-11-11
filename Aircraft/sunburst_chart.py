import matplotlib.pyplot as plt
import  numpy as np

plt.style.use('ggplot')
plt.rcParams['font.sans-serif'] = ['WenQuanYi Micro Hei']
plt.rcParams['legend.fontsize'] = 14

group_names = ['propulsion', 'structures']
group_size = np.array([24, 40]) / 64
subgroup_names = ['engine', 'engine 2', 'fuel ', 'selector', 'test', 'tank', 'pump', 'blade']
subgroup_size = np.array([24, 15, 13, 5, 3, 2, 1, 1])/ 64
labelss = group_names + subgroup_names 
a, b, c = [plt.cm.RdPu, plt.cm.GnBu, plt.cm.Greys]

fig, ax = plt.subplots()
ax.axis('equal') 

mypie2, _ = ax.pie(group_size, radius=1.3-0.3, labels = ['', ''],  colors=[a(0.5), b(0.7)] )
plt.setp(mypie2, width=0.4, edgecolor='white')

mypie, _ = ax.pie(subgroup_size, radius=1.3, labels =  ['','','', '','','', '',''],   colors=[c(0.0), b(0.7), b(0.6), b(0.5), b(0.4), b(0.3), b(0.2), b(0.1)])
plt.setp(mypie, width=0.3, edgecolor='white')

ax.axis('equal')  # Equal aspect ratio ensures a circular pie chart

# Add legend
ax.legend(labelss, loc="best")

# Add legend 
plt.margins(0, 0)

plt.tight_layout() 
plt.show()