# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pylab
import matplotlib as mpl

dir = "output_03_普通"
tStep = 5000

# Z = np.loadtxt(".\\" + dir + "\\task_1_2dVec_P_tStep=" + str(tStep) + ".csv",delimiter=",")
Z = np.loadtxt(".\\task_2_2dVec_P_tStep=" + str(tStep) + ".csv",delimiter=",")
Z = Z.T
for i in xrange(len(Z)):
	for j in xrange(len(Z[0])):
		if Z[i][j] > 0.6:
			Z[i][j] = 0.6
		elif Z[i][j] < -0.6:
			Z[i][j] = -0.6

NX = 401
NY = 201
Tstep2T = 1.0/50.0

x = np.linspace(-10, 30, NX)
y = np.linspace(-10, 10, NY)
print x
print y
X, Y = np.meshgrid(x, y)

# print x
# print y
# exit()

# 一括設定
# http://nanyakan.blogspot.jp/2013/03/matplotlib.html
mpl.rcParams['font.family'] = u'Consolas'
# mpl.rcParams['font.size'] = 14
# params = {'font.family': 'Times New Roman'}
# params = {'font.family': 'Consolas'}
# mpl.rcParams.update(params)
# print mpl.rcParams['font.family']



#plt.contour(X, Y, Z)
# plt.contourf(X, Y, Z)
# plt.contourf(X, Y, Z)
# plt.contour(X, Y, Z, 80)
# colorbarRange = np.linspace(-0.6, 0.6, 30, endpoint=False)
# colorbarRange = np.linspace(-0.6, 0.6, 12*2+1)
# colorbarRange = np.linspace(-0.6, 0.6, 15*2+1)
# colorbarRange = np.linspace(-0.6, 0.6, 12*4+1)
# colorbarRange = np.linspace(-0.6, 0.6, 12*4+1)
colorbarRange = np.linspace(-0.6, 0.6, 12*5+1)
# colorbarRange = np.linspace(-0.3, 0.3, 15*2+1)
print colorbarRange
plt.contour(X, Y, Z, colorbarRange, colors="black",alpha=0.3, linewidths=0.3, linestyles="solid")
plt.contourf(X, Y, Z, colorbarRange, cmap=plt.cm.jet)

# plt.contourf(X, Y, Z, 80, cmap=plt.cm.prism)
colorbarTics = np.linspace(-0.6, 0.6, 6+1)
# pad: プロットエリアとカラーバーの距離
cb = plt.colorbar(label="Pressure", orientation="horizontal", pad=0.1, ticks=colorbarTics)



ax = cb.ax
# font = mpl.font_manager.FontProperties(family='times new roman', style='italic', size=16)
font = mpl.font_manager.FontProperties(family='Consolas')
# text = ax.yaxis.label
text = ax.xaxis.label
# text.set_font_properties(font)


# cb.ax.tick_params(steps=10)
# cb.ax.tick_params(font)
# text = ax.xticks
# text.set_font_properties(font)
# plt.xticks(fontname="Consolas")
# plt.yticks(fontname="Consolas")



rect = pylab.Rectangle((-0.5, -0.5), 1, 1, linewidth=0, facecolor="#777777")
# rect = pylab.Rectangle((-10, -10), 10, 10, linewidth=0, facecolor="#FF0000")
plt.gca().add_patch(rect)

# plt.title("T_Step = " + '{0:05d}'.format(tStep) + ", time = " + str(tStep*Tstep2T))
plt.title("T_Step = " + '{0:05d}'.format(tStep) + ", time = " + '{0:06.1f}'.format(tStep*Tstep2T))
# plt.title("T_Step = " + '{0:05d}'.format(tStep) + ", time = " + '{0:06.1f}'.format(tStep*Tstep2T), fontname="Consolas")

#plt.contour(X, Y, Z, levels=[-0.4, -0.2 ,0, 0.2, 0.4, 0.6, 0.8, 1.0])

plt.gca().set_aspect('equal')
filename = "output_Tstep=" + '{0:05d}'.format(tStep) + ".png"
plt.savefig(filename,dpi=300)
filename = "output_Tstep=" + '{0:05d}'.format(tStep) + ".eps"
plt.savefig(filename)
filename = "output_Tstep=" + '{0:05d}'.format(tStep) + ".pdf"
plt.savefig(filename)
filename = "output_Tstep=" + '{0:05d}'.format(tStep) + ".svg"
plt.savefig(filename)

print np.max(Z), np.min(Z)

#plt.show()
plt.close()

#plt.pause(interval)



