import matplotlib.pyplot as plt

name_list = [ 'our','simfastfor+delta', 'BP32+delta','FOR']
num_list = [ 28.84,18.65,18.65,32.05]
num_list2 = [11.25,1.15,1.45,1.3]
for p1 in range(len(name_list)):
  for p2 in range(len(name_list)):
    if p2>p1:
      if num_list[p1]>num_list[p2]:
        tmp=num_list[p1]
        num_list[p1]=num_list[p2]
        num_list[p2]=tmp
        tmp=num_list2[p1]
        num_list2[p1]=num_list2[p2]
        num_list2[p2]=tmp
        tmp=name_list[p1]
        name_list[p1]=name_list[p2]
        name_list[p2]=tmp
        
x = list(range(len(num_list)))
total_width, n = 0.8, 2
width = total_width / n
plt.bar(x, num_list, width=width,  label='compression rate', fc='b')
for i,j in enumerate(num_list):
  plt.text(i,j+1,'%s' %round(j,1),ha='center')
for i in range(len(x)):
    x[i] += width
plt.bar(x, num_list2, width=width, label='decompress time', tick_label=name_list, fc='g')
for i,j in enumerate(num_list2):
  plt.text(i+width,j+1,'%s' %round(j,1),ha='center')
plt.title('books',fontsize='large',fontweight='bold')
plt.legend()
plt.show()
