import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

# 读取数据
input_data = sys.argv[1]
output_data = sys.argv[2]
column_names = ['genome name', 'genomic location', 'depth coverage']
data = pd.read_csv(input_data, sep='\t', header=None, names=column_names)

# 获取数据
x = data['genomic location']
y = data['depth coverage']

# 计算平均depth coverage、最大深度和最小深度
average_coverage_depth = y.mean()
maximal_depth = y.max()
minimal_depth = y.min()
total_genome_length = x.iloc[-1]

# 创建绘图
plt.figure(figsize=(12, 6))
plt.plot(x, y, color='blue', linewidth=1, linestyle='-')

# 美化图形
plt.title('Sequencing Depth and Coverage Map', fontsize=18, pad=20)
plt.xlabel('Genomic Position (bp)', fontsize=14)
plt.ylabel('Sequencing Depth (×)', fontsize=14)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.grid(True, linestyle='--', alpha=0.5)

# 添加文本
total_genome_length_str = f'{total_genome_length:,} bp'
text_str_1 = f'(1) Total genome length = {total_genome_length_str}'
text_str_3 = f'(3) Maximal depth = {maximal_depth} ×'
text_str_2 = f'(2) Average depth = {average_coverage_depth:.2f} ×'
text_str_4 = f'(4) Minimal depth = {minimal_depth} ×'

plt.gca().annotate(text_str_1, fontsize=12, xy=(0, -0.15), xycoords='axes fraction', xytext=(0, 0), textcoords='offset points', ha='left', va='top')
plt.gca().annotate(text_str_3, fontsize=12, xy=(0, -0.20), xycoords='axes fraction', xytext=(0, 0), textcoords='offset points', ha='left', va='top')
plt.gca().annotate(text_str_2, fontsize=12, xy=(0.5, -0.15), xycoords='axes fraction', xytext=(0, 0), textcoords='offset points', ha='left', va='top')
plt.gca().annotate(text_str_4, fontsize=12, xy=(0.5, -0.20), xycoords='axes fraction', xytext=(0, 0), textcoords='offset points', ha='left', va='top')

# 调整布局以适应文本
plt.subplots_adjust(bottom=0.2)

# 保存图形
plt.savefig(output_data + ".jpg", dpi=600)

# 显示图形
plt.show()

