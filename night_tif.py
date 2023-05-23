from matplotlib import pyplot as plt
import matplotlib
from datetime import timedelta,datetime
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cartopy.io.shapereader as shpreader
import cartopy.crs as ccrs
import gdal,os,shapefile
from gdal import gdalconst

plt.rcParams['font.sans-serif'] = ['STKAITI'] #字体，可修改
def rmfile(fileName):
    ''' 删除文件 '''
    if os.path.exists(fileName):
        os.remove(fileName)
        status = "success"
        message = "ok"
    else:
        status = "warning"
        message = "%s is not exist!" % fileName
    return status, message

def make_dir(path):
    ''' 创建文件夹 '''
    if not os.path.exists(path):
        os.makedirs(path)

def GetTifInfo(tifFile):
    '''获取tif数据属性'''
    dataset = gdal.Open(tifFile, gdal.GA_ReadOnly)  # 打开文件
    #Projection = dataset.GetProjection() # 获得投影信息
    GeoTransform = dataset.GetGeoTransform() # 获得空间参考
    Width = dataset.RasterXSize #列数
    Height = dataset.RasterYSize #行数
    data = dataset.ReadAsArray() # 将数据写成数组，对应栅格矩阵
    resampleFactor = [Width,Height,GeoTransform]
    dataset = None
    return resampleFactor,data

def ProTif(path, fileNames, savePath, outName):
    allTif = []
    for fileName in fileNames:
        dataset = gdal.Open(os.path.join(path,fileName), gdal.GA_ReadOnly)  
        data = dataset.ReadAsArray()
        allTif.append(data)
    sumData = sum(allTif)
    # 创建重采样后的栅格
    make_dir(savePath)
    outFile = os.path.join(savePath, outName)
    
    out_rows = dataset.RasterYSize
    out_columns = dataset.RasterXSize
    num_bands = dataset.RasterCount
    # 创建输出数据集
    gtiff_driver = gdal.GetDriverByName('GTiff')
    out_ds = gtiff_driver.Create(outFile,out_columns, out_rows, num_bands)
    # 像素宽度、像素高度加倍
    out_ds.SetProjection(dataset.GetProjection())
    out_ds.SetGeoTransform(dataset.GetGeoTransform())
    # 将数据写入输出光栅
    out_band = out_ds.GetRasterBand(1)
    out_band.WriteArray(sumData)
    # 将数据写入栅格中
    #out_ds.WriteRaster(0, 0, out_columns, out_rows, sumData)
    out_ds = None


def GdalReprojectImage(tifFile, resampleFactor, savePath, methods):
    '''重采样'''
    dataset = gdal.Open(tifFile, gdal.GA_ReadOnly)
    srcProjection = dataset.GetProjection()
    srcGeoTransform = dataset.GetGeoTransform()
    srcWidth = dataset.RasterXSize
    srcHeight = dataset.RasterYSize
    srcBandCount = dataset.RasterCount
    srcNoDatas = [
		dataset.GetRasterBand(bandIndex).GetNoDataValue()
		for bandIndex in range(1, srcBandCount+1)
	]
    srcBandDataType = dataset.GetRasterBand(1).DataType
    srcFileName = os.path.basename(tifFile)
    name = os.path.splitext(srcFileName)[0]
    # 创建重采样后的栅格
    outFileName = name + ".tif"
    make_dir(savePath)
    outFilePath = os.path.join(savePath, outFileName)
    driver = gdal.GetDriverByName('GTiff')
    geoTransforms = list(srcGeoTransform)
    if isinstance(resampleFactor,list):
        outWidth = resampleFactor[0]
        outHeight = resampleFactor[1]
        geoTransforms[1] = resampleFactor[2][1]
        geoTransforms[5] = resampleFactor[2][5]
    elif isinstance(resampleFactor,int):
        outWidth = int(srcWidth * resampleFactor)
        outHeight = int(srcHeight * resampleFactor)
        geoTransforms[1] = geoTransforms[1]/resampleFactor
        geoTransforms[5] = geoTransforms[5]/resampleFactor
    else:
        print('重采样因子resampleFactor 错误')
    outDataset = driver.Create(
        outFilePath,#输出栅格
		outWidth,
		outHeight,
		srcBandCount,
		srcBandDataType
	)
    outGeoTransform = tuple(geoTransforms)
    outDataset.SetGeoTransform(outGeoTransform)
    outDataset.SetProjection(srcProjection)
    for bandIndex in range(1, srcBandCount+1):
        band = outDataset.GetRasterBand(bandIndex)
        band.SetNoDataValue(srcNoDatas[bandIndex-1])
    gdal.ReprojectImage(
		dataset,
		outDataset,
		srcProjection,
		srcProjection,
		methods,
	)
    
    resampleData = outDataset.ReadAsArray()
    del outDataset
    return resampleData

def TiffCutByShp(tifFile, shpfile, savePath):
    dataset = gdal.Open(tifFile, gdal.GA_ReadOnly)
    datashp = shapefile.Reader(shpfile)
    srcFileName = os.path.basename(tifFile)
    name = os.path.splitext(srcFileName)[0]
    outFileName = name + ".tif"
    outFilePath = os.path.join(savePath, outFileName)
    outDataset = gdal.Warp(outFilePath,
                           dataset,
                           format='GTiff', 
                           cutlineDSName = shpfile,
                           outputBounds = datashp.bbox,
                           dstNodata = 0)  
    cutData = outDataset.ReadAsArray()
    cutWidth = outDataset.RasterXSize
    cutHeight = outDataset.RasterYSize
    GeoTransform = outDataset.GetGeoTransform()
    resampleFactor = [cutWidth, cutHeight, GeoTransform]
    outDataset = None
    return resampleFactor,cutData

def landclass(data, grades):
    try:
        data = np.where(np.isnan(data), 0, data)
        data = np.where(data == -9999, 0, data)
    except:
        pass
    for one in grades:
        data = np.where(data == int(one), grades[one], data)
    return data
"""
def plot_show(tifdata, sets, shpfile, pngName):
    #fig = plt.figure()                
    #ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    sf = shapefile.Reader(shpfile)
    border_shape = sf
    border = border_shape.shapes()
    for i in range(len(border)):
        border_points = border[i].points
        x, y = zip(*border_points)
        plt.plot(x, y,'k-')  # x横坐标 y纵坐标 ‘k-’线性为黑色
    lonmin, lonmax = sf.bbox[0], sf.bbox[2]
    latmin, latmax = sf.bbox[1], sf.bbox[3]
    lon_xticks = np.array(['%.2f'%x for x in np.linspace(lonmin,lonmax,5)],dtype=np.float64)
    lat_yticks = np.array(['%.2f'%x for x in np.linspace(latmin,latmax,5)],dtype=np.float64)
    norm = matplotlib.colors.Normalize(vmin=sets['bound'][0],vmax=sets['bound'][1])  # 设置colorbar显示的最大最小值
    im = plt.imshow(tifdata, 
                    cmap=sets['cmap'],
                    alpha=sets['alpha'],
                    norm=norm,
                    extent = (lonmin, lonmax, latmin, latmax),
                    )  # 绘图，配色
    plt.title(pngName.split('_')[0] + sets['title'])
    plt.xticks(lon_xticks,[str(x)+'°E' for x in list(lon_xticks)],rotation=0,fontsize=10,fontweight='light')
    plt.yticks(lat_yticks,[str(y)+'°N' for y in list(lat_yticks)],fontsize=10,fontweight='light')
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", size="5%", pad="3%")
    cbar = plt.colorbar(im, cax=cax)#图例
    cbar.set_label('kg/h',fontsize=12),#loc='center  
    shpplot = shpreader.Reader(shpfile)
    im.add_geometries(shpplot.geometries(), crs=ccrs.PlateCarree(), linewidths=0.5, edgecolor='k',facecolor='none')
    make_dir(sets['pngPath'])
    plt.savefig(os.path.join(sets['pngPath'],pngName), dpi=300)
    plt.close()
"""  
def plot_show(tifdata, sets, shpfile, pngName):
    fig = plt.figure()                
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    sf = shapefile.Reader(shpfile)
    lonmin, lonmax = sf.bbox[0], sf.bbox[2]
    latmin, latmax = sf.bbox[1], sf.bbox[3]
    lon_xticks = np.array(['%.2f'%x for x in np.linspace(lonmin,lonmax,sets['xticknum'])],dtype=np.float64)
    lat_yticks = np.array(['%.2f'%x for x in np.linspace(latmin,latmax,sets['yticknum'])],dtype=np.float64)
    im = ax.imshow(tifdata, 
                   cmap=sets['cmap'], # 配色
                   alpha=sets['alpha'], #透明度
                   vmin=sets['bound'][0],vmax=sets['bound'][1],# 设置colorbar显示的最大最小值
                   extent=(lonmin, lonmax, latmin, latmax),
                   origin='upper',
                   )  
    cbar = plt.colorbar(im, ax=ax, shrink=sets['shrink'], 
                        pad=sets['pad'], orientation=sets['orientation'],
                        extend="both",)
    cbar.ax.tick_params(labelsize='small')
    cbar.set_label(label=sets['unit'],fontsize=sets['unitsize'],
                   fontname='Times New Roman', fontweight='bold')
    
    shpplot = shpreader.Reader(shpfile)
    ax.add_geometries(shpplot.geometries(), crs=ccrs.PlateCarree(), 
                      linewidths=0.5, edgecolor='k',facecolor='none') 
    plt.xticks(lon_xticks,[str(x)+'°E' for x in lon_xticks.tolist()],
               rotation=sets['rotation'],fontsize=sets['xyfontsize'],
               fontname='Times New Roman', fontweight='bold')
    plt.yticks(lat_yticks,[str(y)+'°N' for y in lat_yticks.tolist()],
               fontsize=sets['xyfontsize'],fontname='Times New Roman', 
               fontweight='bold')
    plt.title(pngName.split('_')[0] + sets['title'], fontsize=sets['titlesize'], 
              fontname='Times New Roman', weight='bold')
    make_dir(sets['pngPath'])
    plt.savefig(os.path.join(sets['pngPath'],pngName), dpi=sets['dpi'])
    plt.close(fig)    
    

if __name__ == '__main__':
    """===========================================修改参数==========================================="""
    begTime = '20190112'  #起始时间
    ndays = 1  #时长：day
    grades = {'51':5, '52':6,'53':7} #土地数据分类等级,并赋值
    method = gdalconst.GRA_Bilinear
    #gdalconst.GRA_NearestNeighbour：最邻近插值
    #gdalconst.GRA_Bilinear:双线性插值
    #gdalconst.GRA_Cubic:三次插值
    #gdalconst.GRA_CubicSpline:三次样条插值
    #gdalconst.GRA_Lanczos:lanczos
    #gdalconst.GRA_Average:average
    #gdalconst.GRA_Mode:mode
    
    ###输入数据路径###
    landFile = 'D:/co2_ghy/ld2020_1000m_gs84.tif' #土地利用数据，wgs84
    shpfile = 'D:/co2_ghy/shpfile/province.shp' #文件
    lightPath = 'D:/co2_ghy' #夜光数据路径
    ###输出数据路径（不用新建文件夹,会自动识别并创建）###
    TiffCutPath = 'D:/co2_ghy/TiffCut' #tif裁剪保存路径
    resamplePath = 'D:/co2_ghy/reshape' #重采样tif存放路径
    ###画图参数配置###
    sets = dict(title = '__CO2 emission inventory of 100m', #图片标题
                titlesize = 12,#图片标题字体大小
                cmap = 'bwr', #颜色：'bwr','rainbow','viridis','Blues','YlGnBu','Greens'
                pngPath = 'D:/co2_ghy/fig', #图片存放路径
                alpha = 0.8, #透明度
                bound = [0,50],# 设置colorbar显示的最大最小值
                unit = 'kg/h',  #色带标记（例如单位）
                unitsize = 12, #色带标记字体大小
                shrink = 0.64, #缩放参数shrink，从0-1，整个色条将会按照输入值被缩放0.6
                pad = 0.1, #色条与子图的间距,0.01,0.1
                orientation = 'horizontal', #控制色条时横纵方向，"vertical",'horizontal'
                xyfontsize = 10,#xy轴刻度字体大小
                xticknum = 5, #x轴刻度显示个数
                rotation = 0, #x轴刻度旋转度数
                yticknum = 5, #y轴刻度显示个数
                dpi = 300,#图片分辨率
                )
    """============================================================================================="""
    begTime = datetime.strptime(begTime, "%Y%m%d")
    make_dir(TiffCutPath)
    resampleFactor, cutData = TiffCutByShp(landFile, shpfile, TiffCutPath)#获取tif文件裁剪后的分辨率
    #resampleFactor,tifData = GetTifInfo(landFile)#获取tif文件的分辨率
    #resampleFactor = 2
    classData = landclass(cutData, grades)
    for i in range(ndays):
        thisTime = begTime + timedelta(days=i)
        fileName = thisTime.strftime('%Y-%m-%d_notnight.tif')
        if os.path.exists(os.path.join(lightPath,fileName)):
            other, light_cutData = TiffCutByShp(os.path.join(lightPath,fileName), shpfile, TiffCutPath)
            resampleData = GdalReprojectImage(os.path.join(TiffCutPath,fileName), resampleFactor, resamplePath, method)
            emisData = resampleData * classData
            pngName = os.path.splitext(fileName)[0] + ".png"
            plot_show(emisData, sets, shpfile, pngName)
        else:
            continue