import matplotlib.pyplot as plt

# 创建数据
x = [5, 10, 20, 40, 80, 160, 320, 640]  # x轴数据
y = [0.75, 0.76, 0.77, 0.78, 0.81, 0.86, 0.95, 1.13]  # y轴数据

# 绘制折线图
plt.plot(x, y)

# 添加标题和坐标轴标签
plt.title("折线图示例")
plt.xlabel("X轴")
plt.ylabel("Y轴")

# 显示图形
plt.show()
