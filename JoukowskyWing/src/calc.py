# coding: UTF-8

"""
文字コードはutf-8N
改行コードは CRCF
Windows cmdではユニコード文字列を使うこと．
print u"にほんご"
print "にほんご"

2016/06/28

【コーディング規約】
ユニコード文字列は""
普通の文字列は''
で書く

ハンガリー記法
c : 複素数
f : float
i : int
"""


import sys
import os
import numpy as lNP
import pylab as lPL
import cmath as lcM
import math as lfM

DEBUG = 0
# show or save
IS_SHOW_GRAPH = 0
fPI = lfM.pi

def main():
    fA = 0.60
    fV0 = 1
    fALPHA = deg2rad(5.0)
    fBETA = deg2rad(20.0)
    fC = 0.5
    # fGAMMA = 5
    fGAMMA = 4 * fPI * fV0 * fA * lfM.sin(fALPHA + fBETA)

    # 極座標のメッシュ分割数
    iR_GRID = 201
    iTHETA_GRID = 3601   #  一周回すために，360 + 1
    # iTHETA_GRID = 361   #  一周回すために，360 + 1
    fR_MIN = fA         # r方向は，翼型の縁からスタート
    fR_MAX = 3

    ###########################################
    # 翼型の変換先cZetaを求める
    fAirfoilXis = [0] * iTHETA_GRID            # 翼型
    fAirfoilEtas = [0] * iTHETA_GRID           # 翼型
    for iTheta in range(iTHETA_GRID):
        cLargeZ = fA * lcM.exp( (deg2rad(iTheta) + fALPHA) * 1.0j ) + fA * lcM.exp( (fPI - fBETA) * 1.0j ) + fC
        cZeta = cLargeZ + fC**2 / cLargeZ
        fAirfoilXis[iTheta] = cZeta.real
        fAirfoilEtas[iTheta] = cZeta.imag

    ###########################################
    # 流線ψの変換先cZetaを求める
    # そのためにまずZ=X+iYグリッドをもとめる．
    # 等高線
    # http://d.hatena.ne.jp/y_n_c/20091122/1258904025
    """
    # data format
    x =                     y =                             f(z) =
    [                       [                               [
    [1, 1, 1, 1, 1]         [0.1, 0.2, 0.3, 0.4, 0.5]       [0.1, 0.2, 0.3, 0.1, 0.3]
    [2, 2, 2, 2, 2]         [0.1, 0.2, 0.3, 0.4, 0.5]       [0.4, 0.5, 0.1, 0.9, 0.5]
    [3, 3, 3, 3, 3]         [0.1, 0.2, 0.3, 0.4, 0.5]       [0.5, 0.2, 0.2, 0.4, 0.1]
    [4, 4, 4, 4, 4]         [0.1, 0.2, 0.3, 0.4, 0.5]       [0.6, 0.3, 0.3, 0.4, 0.1]
    [5, 5, 5, 5, 5]         [0.1, 0.2, 0.3, 0.4, 0.5]       [0.1, 0.2, 0.7, 0.4, 0.5]
    ]                       ]                               ]

    → r
    ↓ Θ
    """
    # 変数初期化
    # fXs = [[0 for col in range(iR_GRID)] for row in range(iTHETA_GRID)]
    # fYs = [[0 for col in range(iR_GRID)] for row in range(iTHETA_GRID)]
    fXis = [[0 for col in range(iR_GRID)] for row in range(iTHETA_GRID)]
    fEtas = [[0 for col in range(iR_GRID)] for row in range(iTHETA_GRID)]
    fPsis = [[0 for col in range(iR_GRID)] for row in range(iTHETA_GRID)]
    fCps = [[0 for col in range(iR_GRID)] for row in range(iTHETA_GRID)]
    fCpWallXis = []             # 壁面のξ座標
    fCpWalls = []               # 壁面のξ座標に対するCp
    # 圧力の積分用
    fCpXiSum = 0.0
    fCpEtaSum = 0.0
    fPreXi = 0.0
    fPreEta = 0.0
    fCpAvg = 0.0
    fPreCp = 0.0
    # 風圧中心用
    fXcpBunbo = 0.0
    fXcpBunshi = 0.0
    # 空力中心用
    fXisForfCms = [i*0.01 for i in range(-55,-34)]        # [-1.00, -0.95, -0.90, ... , 0.00]
    fCms = [0.0] * len(fXisForfCms)


    dTheta = 360.0 / (iTHETA_GRID - 1)
    dR = 1.0*(fR_MAX - fR_MIN) / (iR_GRID - 1)

    for iRow in range(iTHETA_GRID):
        fTheta = deg2rad(0.0 + dTheta * iRow)
        if (DEBUG == 1): print rad2deg(fTheta)
        for iCol in range(iR_GRID):
            # xyグリッド
            fR = fR_MIN + dR * iCol
            fX = fR * lfM.cos(fTheta)
            fY = fR * lfM.sin(fTheta)
            cZ = complex(fX, fY)
            # Z計算
            cLargeZ = cZ * lcM.exp(fALPHA * 1.0j) + fA * lcM.exp( (fPI - fBETA) * 1.0j ) + fC

            # f,ψを求める
            # cF = fV0 * (cLargeZ + fA**2 / cLargeZ) + ((fGAMMA * 1.0j) / (2*fPI)) * lcM.log(cZ)
            cF = fV0 * (cZ + fA**2 / cZ) + ((fGAMMA * 1.0j) / (2*fPI)) * lcM.log(cZ)
            cZeta = cLargeZ + fC**2 / cLargeZ
            fXi = cZeta.real
            fEta = cZeta.imag
            fXis[iRow][iCol] = fXi
            fEtas[iRow][iCol] = fEta
            fPsis[iRow][iCol] = cF.imag

            # Cpを求める
            cdFdZeta = ( fV0 * lcM.exp(-fALPHA * 1.0j) * (1 - (fA**2 / cZ**2) + ((fGAMMA * 1.0j) / (2*fPI*cZ))) ) / (1 - (fC**2 / cLargeZ**2))
            fCp = 1.0 - ( (cdFdZeta.real)**2 + (cdFdZeta.imag)**2 ) / fV0**2
            fCps[iRow][iCol] = fCp

            # 壁面Cpを抽出する
            if (fR == fR_MIN):
                fCpWallXis.append(fXi)
                fCpWalls.append(fCp)
                # 壁面Cpの積分
                if (iRow == 0):
                    pass
                else:
                    fdXiWall = fXi - fPreXi
                    fdEtaWall = fEta - fPreEta
                    fCpAvg = (fCp + fPreCp) / 2.0
                    fCpXiSum += - fCpAvg * fdEtaWall
                    fCpEtaSum += - fCpAvg * (-fdXiWall)
                    # print fCpXiSum, fCpEtaSum
                    # 風圧中心計算
                    fXiWallAvg = (fXi + fPreXi)  /2.0
                    fEtaWallAvg = (fEta + fPreEta)  /2.0
                    fXcpBunshi += fXiWallAvg * ( - fCpAvg * (-fdXiWall) ) - fEtaWallAvg * ( - fCpAvg * fdEtaWall )
                    fXcpBunbo += - fCpAvg * (-fdXiWall)
                    # 空力中心計算
                    for i in range(len(fXisForfCms)):
                        fCms[i] += ( fXisForfCms[i] - fXiWallAvg ) * ( - fCpAvg * (-fdXiWall) ) + fEtaWallAvg * ( - fCpAvg * fdEtaWall )
                fPreXi = fXi
                fPreEta = fEta
                fPreCp = fCp


    # 日本語フォント
    # http://qiita.com/canard0328/items/a859bffc9c9e11368f37
    from matplotlib.font_manager import FontProperties
    # fp = FontProperties(fname=r'C:\WINDOWS\Fonts\YuGothic.ttf', size=14)
    fp = FontProperties(fname=r'C:\WINDOWS\Fonts\YuGothR.ttc', size=14)

    lPL.figure(figsize=(9, 6), dpi=80)
    # lPL.figure(dpi = 800)
    # lPL.axis([-2.0, 2.0, -1.0, 2.5])
    lPL.axis([-2.0, 2.0, -1.0, 1.5])
    lPL.gca().set_aspect('equal', adjustable='box')    # 縦横比を揃える
    lPL.plot(fAirfoilXis, fAirfoilEtas, "black", linewidth=2)
    # lPL.plot(x2,y2, "blue")
    # lPL.contour(fXis, fEtas, fPsis, 80)    # 等高線
    lPL.contour(fXis, fEtas, fPsis, [i*0.05 for i in range(-19,37)])    # 等高線
    # lPL.contour(fXis, fEtas, fPsis, [i*0.05 for i in range(-21,35)])    # 等高線 fALPHA = -5 の時用
    # lPL.colorbar(label='psi', family='YuGothR')          # カラーバー
    lPL.colorbar(label='psi')          # カラーバー
    # lPL.colorbar()          # カラーバー
    # lPL.label('psi')
    lPL.spring()            # 連続変化色
    # lPL.xlabel(u'ξ', fontproperties=fp)
    # lPL.ylabel(u'η', fontproperties=fp)
    lPL.xlabel('xi')
    lPL.ylabel('eta')
    if (IS_SHOW_GRAPH == 1):
        lPL.show()
    else:
        lPL.savefig('psi.png')
        lPL.savefig('psi.svg')
        lPL.savefig('psi.eps')
        lPL.savefig('psi.ps')
        lPL.savefig('psi.pdf')
    # return

    lPL.figure(figsize=(9, 6), dpi=80)
    # lPL.figure(dpi = 800)
    # lPL.axis([-2.0, 2.0, -1.0, 2.5])
    # lPL.axis([-2.0, 2.0, -1.0, 1.5])
    lPL.gca().set_aspect('equal', adjustable='box')    # 縦横比を揃える
    lPL.plot(fAirfoilXis, fAirfoilEtas, "black", linewidth=2)
    # lPL.plot(x2,y2, "blue")
    lPL.contour(fXis, fEtas, fCps, 80)    # 等高線
    # lPL.contour(fXis, fEtas, fCps, [i*0.05 for i in range(-19,37)])    # 等高線
    lPL.colorbar(label='Cp')          # カラーバー
    lPL.spring()            # 連続変化色
    # lPL.xlabel(u'ξ', fontproperties=fp)
    # lPL.ylabel(u'η', fontproperties=fp)
    lPL.xlabel('xi')
    lPL.ylabel('eta')
    if (IS_SHOW_GRAPH == 1):
        lPL.show()
    else:
        lPL.savefig('Cp.png')
        lPL.savefig('Cp.svg')

    lPL.figure(figsize=(9, 6), dpi=80)
    # lPL.figure(dpi = 800)
    # lPL.axis([-2.0, 2.0, -1.0, 2.5])
    lPL.axis([-1.5, 1.5, -3.5, 1.5])
    lPL.gca().set_aspect('equal', adjustable='box')    # 縦横比を揃える
    lPL.plot(fAirfoilXis, fAirfoilEtas, "black", linewidth=2)
    # lPL.plot(x2,y2, "blue")
    lPL.scatter(fCpWallXis, fCpWalls, s=10, marker='.', linewidths=0)
    lPL.grid(True)
    # lPL.xlabel(u'ξ', fontproperties=fp)
    lPL.xlabel('xi')
    lPL.ylabel('Cp')
    if (IS_SHOW_GRAPH == 1):
        lPL.show()
    else:
        lPL.savefig('Cp_wall.png')
        lPL.savefig('Cp_wall.svg')

    fCpXiSum = fCpXiSum / (4.0 * fC)
    fCpEtaSum = fCpEtaSum / (4.0 * fC)
    fCl = fCpEtaSum * lfM.cos(fALPHA) - fCpXiSum * lfM.sin(fALPHA)
    fCd = fCpXiSum * lfM.cos(fALPHA) + fCpEtaSum * lfM.sin(fALPHA)

    ps("\n")
    ps("CL = ")
    pr(fCl)
    ps("CD = ")
    pr(fCd)

    fXcp = fXcpBunshi / fXcpBunbo
    ps("Xcp = ")
    pr(fXcp)

    # ρ= 1 とした．
    for i in range(len(fCms)):
        fCms[i] = fCms[i] / (0.5 * 1.0 * fV0**2 * (4 * fC)**2)
    outputFile = open('Xcp-Cm.csv', "a")
    outputFile.write('%.6f' % rad2deg(fALPHA))
    for i in range(len(fCms)):
        outputFile.write(',')
        outputFile.write('%.6f' % fCms[i])
    outputFile.write('\n')
    outputFile.close()


    return

    """
    from pylab import *
    from numpy import *

    x = arange(0, 10.1, 0.1)
    y = arange(0, 10.1, 0.1)
    x = arange(0, 10.1, 1)
    y = arange(0, 10.1, 1)
    # x
    # y

    X, Y = meshgrid(x, y)

    Z = cos(X) + cos(Y)
    Z

    contour(X, Y, Z)
    colorbar()
    spring()

    show()
    """



def deg2rad(fDeg):
    return fDeg * fPI / 180

def rad2deg(fRad):
    return fRad * 180 / fPI

def ps(output):
    sys.stdout.write(str(output))
    sys.stdout.flush()

# ユニコード文字列はstr()にいれるとバグったので
def pu(output):
    sys.stdout.write(output)
    sys.stdout.flush()

def pr(string):
    print (string)
    sys.stdout.flush()

def GetStringTime():
    import datetime # datetimeモジュールのインポート
    import locale

    d = datetime.datetime.today()   # モジュール名.クラス名.メソッド名
    stringTime = str(d.year)
    stringTime += '.'
    stringTime += '{0:02d}'.format(d.month)
    stringTime += '.'
    stringTime += '{0:02d}'.format(d.day)
    stringTime += '-'
    stringTime += '{0:02d}'.format(d.hour)
    stringTime += '.'
    stringTime += '{0:02d}'.format(d.minute)
    stringTime += '.'
    stringTime += '{0:02d}'.format(d.second)
#    stringTime += '.'
#    stringTime += '{0:02d}'.format(d.microsecond / 10000)
    return stringTime

if __name__ == '__main__':
    main()




