import matplotlib.pyplot as plt
from scipy import interpolate

def interpolate_points(x_points, y_points):
    # 线性插值
    interpolation_function = interpolate.interp1d(x_points, y_points, kind='linear')

    # 返回一个函数来计算任意点的函数值
    def calculate_function_value(x):
        return interpolation_function(x)

    return calculate_function_value

def plot_line_graph(x, y1, y2, y3, y4, x_label, y_labels, save_file, text_dots = []):
    plt.figure(figsize=(10, 4))  # 设置图形尺寸

    plt.plot(x, y1, label=y_labels[0])
    plt.plot(x, y2, label=y_labels[1])
    plt.plot(x, y3, label=y_labels[2])
    plt.plot(x, y4, label=y_labels[3])
    
    for dot in text_dots:
        plt.text(25000, dot+10, '%.1f' % (dot))
    # plt.text(25000, 260, '220')  # 标记数据点
    # plt.text(25000, 880, '840')  # 标记数据点

    plt.xlabel(x_label)  # 自变量标签
    plt.ylabel('Time (ms)')  # 因变量的通用标签

    plt.legend()  # 显示图示

    # plt.title('折线图')  # 图表标题
    plt.grid(True)  # 显示网格
    plt.savefig(save_file)
    
    plt.show()

# 示例数据
x = [2340, 4680, 9360, 18720, 37440]  # 自变量数据
y1 = [287.5, 304.6, 320.6, 433.5, 587.8]  # 因变量1
y2 = [142.7, 150.9, 168.6, 204.2, 280.3]  # 因变量2
y3 = [27, 45.2, 94.3, 153.2, 325.1]  # 因变量3
y4 = [13.6, 24.9, 45.1, 87.7, 169.2]  # 因变量4
for i in y1:
    print(i)
for i in y2:
    print(i)
x_label = 'original circuit size |C|'  # 自变量标签
y_labels = ['PI_ILC Prove', 'PI_ILC Verify', 'CP_ILC Prove', 'CP_ILC Verify']  # 因变量标签

f_y1_value = interpolate_points(x, y1)(25000)
f_y2_value = interpolate_points(x, y2)(25000)
f_y3_value = interpolate_points(x, y3)(25000)
f_y4_value = interpolate_points(x, y4)(25000)
text_dots = [f_y1_value, f_y2_value, f_y3_value, f_y4_value]
print(text_dots)

# 调用函数绘制折线图
plot_line_graph(x, y1, y2, y3, y4, x_label, y_labels, "case1.png", text_dots)

plt.close()

# 示例数据
x = [1, 2, 4, 8, 16, 32]  # 自变量数据
y1 = [446, 765.2, 1286, 2420, 4808, 10509]  # 因变量1
y2 = [216, 365, 643.9, 1206, 2371, 4891]  # 因变量2
y3 = [155.7, 158.7, 164.2, 159.8, 170.8, 167.9]  # 因变量3
y4 = [88.9, 86.5, 89.2, 89.3, 97.5, 93.6]  # 因变量4
x_label = 'commitment num'  # 自变量标签
y_labels = ['PI_ILC Prove', 'PI_ILC Verify', 'CP_ILC Prove', 'CP_ILC Verify']  # 因变量标签

# 调用函数绘制折线图
plot_line_graph(x, y1, y2, y3, y4, x_label, y_labels, "case2.png")
plt.close()
