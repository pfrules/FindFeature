# function PixelPoints=subpixelDetection(rawImageFolder)
        基于杨洋所写使用形态学和背景相减，加上自己的函数:
## result=inmyway(img,centers, radii);
        检测到的圆心半径为粗略估计，然后以拟合的方式选择梯度最大的亚像素或者灰度值相等的亚像素
### 灰度值相等
        实验证明灰度值相等的方法对于实际图像会出现缺少圆的边缘缺少的情况，因为钢球图像在两边出现穿透效果不好的情况。很难找到合适的灰度值，从而使得部分检测到的边缘缺少一部分。
### 梯度最大
        使用梯度最大话，由于图像本身的原因，边缘不够明显，对于灰度值拟合的方法加上梯度最大话效果并不好，很多情况下出现检测特征点跑飞的情况。效果更次于灰度值相等的方法


# PixelPoints=newSubpixelDetection(rawImageFolder)
    另一种检测方法，先用canny算子估算边缘，然后做其他操作
#### 注意：
        由于干扰项较多，对于canny检测到的边缘，即使使用霍夫变换也很难提取出有用的边缘。出现较多干扰。为此需要结合形态学处理进行粗略定位圆的大致位置。
 找到边缘后，排除圆外特征点，方法如下：
 ### 一、
 在园内任选一点，在360度方向上计算点到特征点距离，建立特殊模型，拟合模型曲线，递归调用删除模型外点。
 ### 二、
 直接进行圆拟合，递归调用逐渐删除距离较远的点。

 # PixelPoints=otherSubpixelDetection(rawImageFolder)
 同样使用边缘检测，但是，在圆检测过程中速度较慢。形态学处理后删掉圆检测部分直接进行最大连通域检测