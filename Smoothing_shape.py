#!/usr/bin/env python
# _*_ coding: utf-8 _*_
__author__ = 'zyx'
__date__ = '2022/2/7 18:06'
import math
import sys
import fiona
import osgeo
import shapefile
from shapely.geometry import shape,Polygon,mapping,LinearRing
import numpy as np
def smoothing_shape(infile,shapely,smoothing_order,cell_size,outfile):
    cs_multi=[1.0, 0.7071, 0.6533, 0.6407, 0.6376, 0.6369, 0.6367, 0.6366]
    if smoothing_order>len(cs_multi):
        print("the smoothing order must smaller than 5")
        return
    for i in range(len(infile)):
        geom=shape(infile[i].get('geometry'))
        attributes=infile[i].get('properties')
        a=attributes.get('RGIId')
        if len(geom.interiors)==0:
            point=geom.exterior
            x=[]
            y=[]
            for j in range(0,len(point.xy[0])):
                # x.append(round(point.xy[0][j],4))
                # y.append(round(point.xy[1][j],4))
                x=list(point.xy[0])
                y=list(point.xy[1])
            n_parts=len(shapely.shape(i).parts)
            parts=shapely.shape(i).parts
            bSmmothed=0
            for j in range(0,smoothing_order):
                cs=float(cell_size)*cs_multi[j]
                result,newparts,newx,newy=Smoothing_line(x,y,n_parts,parts,cs)
                if result>0:
                    break
                else:
                    bSmmothed=1
                    x=newx
                    y=newy
                    parts=newparts
            line = []
            if bSmmothed:
                for k in range(len(x)):
                    line.append([x[k], y[k]])
                polygon = Polygon(line)
                polygon = mapping(polygon)
                feature = {
                    'geometry': polygon,
                    'properties': {
                        'RGIId':a,'GLIMSId':attributes.get('GLIMSId'),
                        'BgnDate':attributes.get('BgnDate'),
                        'EndDate':attributes.get('EndDate'),
                        'CenLon':attributes.get('CenLon'),
                        'CenLat':attributes.get('CenLat'),
                        'O1Region':attributes.get('O1Region'),
                        'O2Region':attributes.get('O2Region'),
                        'Area':attributes.get('Area'),
                        'Zmin':attributes.get('Zmin'),
                        'Zmax':attributes.get('Zmax'),
                        'Zmed':attributes.get('Zmed'),
                        'Slope':attributes.get('Slope'),
                        'Aspect':attributes.get('Aspect'),
                        'Lmax':attributes.get('Lmax'),
                        'Status':attributes.get('Status'),
                        'Connect':attributes.get('Connect'),
                        'Form':attributes.get('Form'),
                        'TermType':attributes.get('TermType'),
                        'Surging':attributes.get('Surging'),
                        'Linkages':attributes.get('Linkages'),
                        'Name':attributes.get('Name')
                    },
                }
                outfile.write(feature)
        else:
            x=[]
            y=[]
            exterior_point = geom.exterior
            for k in range(0,len(exterior_point.xy[0])):
                x.append(exterior_point.xy[0][k])
                y.append(exterior_point.xy[1][k])
            # for k in range(len(geom.interiors)):
            #     interiors_point=geom.interiors[k]
            #     for l in range(len(geom.interiors[k].xy[0])):
            #         x.append(interiors_point.xy[0][l])
            #         y.append(interiors_point.xy[0][l])
            n_parts = len(shapely.shape(i).parts)
            parts = list(shapely.shape(i).parts)
            bSmmothed = 0
            for j in range(smoothing_order):
                cs = float(cell_size) * cs_multi[j]
                result, newparts, newx, newy = Smoothing_line(x, y, 1, [0], cs)
                if result > 0:
                    break
                else:
                    bSmmothed = 1
                    x = newx
                    y = newy
                    parts = newparts
            line = []
            if bSmmothed:
                for k in range(len(x)):
                    line.append([x[k], y[k]])
            interiors_line=[]
            interiors_lines=[]
            for k in range(len(geom.interiors)):
                interiors_point=geom.interiors[k]
                x=[]
                y=[]
                for l in range(len(geom.interiors[k].xy[0])):
                    x.append(interiors_point.xy[0][l])
                    y.append(interiors_point.xy[1][l])

                for j in range(smoothing_order):
                    cs = float(cell_size) * cs_multi[j]
                    result, newparts, newx, newy = Smoothing_line(x, y, 1, [0], cs)
                    if result > 0:
                        break
                    else:
                        bSmmothed = 1
                        x = newx
                        y = newy
                        parts = newparts
                if bSmmothed:
                    for m in range(0,len(x)):
                        if len(x)<3:
                            continue
                        interiors_line.append([x[m],y[m]])
                if len(interiors_line)!=0:
                    interiors_lines.append(interiors_line)
                interiors_line=[]
            if bSmmothed:
                try:
                    polygon = Polygon(shell=line,holes=interiors_lines)
                except:
                    print('error')
                polygon = mapping(polygon)
                feature = {
                    'geometry': polygon,
                    'properties': {
                        'RGIId': a, 'GLIMSId': attributes.get('GLIMSId'),
                        'BgnDate': attributes.get('BgnDate'),
                        'EndDate': attributes.get('EndDate'),
                        'CenLon': attributes.get('CenLon'),
                        'CenLat': attributes.get('CenLat'),
                        'O1Region': attributes.get('O1Region'),
                        'O2Region': attributes.get('O2Region'),
                        'Area': attributes.get('Area'),
                        'Zmin': attributes.get('Zmin'),
                        'Zmax': attributes.get('Zmax'),
                        'Zmed': attributes.get('Zmed'),
                        'Slope': attributes.get('Slope'),
                        'Aspect': attributes.get('Aspect'),
                        'Lmax': attributes.get('Lmax'),
                        'Status': attributes.get('Status'),
                        'Connect': attributes.get('Connect'),
                        'Form': attributes.get('Form'),
                        'TermType': attributes.get('TermType'),
                        'Surging': attributes.get('Surging'),
                        'Linkages': attributes.get('Linkages'),
                        'Name': attributes.get('Name')
                    },
                }
                outfile.write(feature)
        if bSmmothed :
            print(a)
def Smoothing_line(x,y,n_parts,parts,cellsize):
    sx_old=[]
    sy_old=[]
    npts=len(x)
    if npts!=len(y):
        return 2
    if npts<=1:
        return 3
    if npts==2:
        sx=[x[0],x[1]]
        sy=[y[0],y[1]]
        newparts=parts
        return 0,newparts,sx,sy
    slp_tol=pow(0.1,5)
    newparts1=parts
    bFrist=1
    for i in range(0,n_parts):
        if bFrist:
            newparts1[i]=0
        else:
            newparts1[i]=len(sx_old)
        length=math.sqrt(math.pow(x[parts[i]]-x[parts[i]+1],2)+math.pow(y[parts[i]]-y[parts[i]+1],2))
        if length-cellsize>slp_tol:
            x01 = (length * x[parts[i]] - cellsize / 2.0 * (x[parts[i]] - x[parts[i] + 1])) / length
            y01 = (length * y[parts[i]] - cellsize / 2.0 * (y[parts[i]] - y[parts[i] + 1])) / length
            x02 = (length * x[parts[i]] - (length- cellsize / 2.0) * (x[parts[i]] - x[parts[i] + 1])) / length
            y02 = (length * y[parts[i]] - (length - cellsize / 2.0) * (y[parts[i]] - y[parts[i] + 1])) / length
            if bFrist:
                sx_old.append(x01)
                sx_old.append(x02)
                sy_old.append(y01)
                sy_old.append(y02)
                # sx_old = [x01, x02]
                # sy_old = [y01, y02]
            else:
                sx_old.append(x01)
                sx_old.append(x02)
                sy_old.append(y01)
                sy_old.append(y02)
                # sx_old = [sx_old,x01, x02]
                # sy_old = [sy_old,y01, y02]
        else:
            if bFrist:
                sx_old.append(float((x[parts[i]] + x[parts[i] + 1]) / 2.0))
                sy_old.append(float((y[parts[i]] + y[parts[i] + 1]) / 2.0))
            else:
                # sx_old = [sx_old,float((x[parts[i]] + x[parts[i] + 1]) / 2.0)]
                # sy_old = [sy_old,float((y[parts[i]] + y[parts[i] + 1]) / 2.0)]
                sx_old.append(float((x[parts[i]] + x[parts[i] + 1]) / 2.0))
                sy_old.append(float((y[parts[i]] + y[parts[i] + 1]) / 2.0))
        bFrist=0

        if i== n_parts-1:
            limit=npts-2
        else:
            limit=parts[i+1]-2
        for j in range(parts[i]+1,limit+1):
            length=math.sqrt(math.pow(x[j]-x[j+1],2)+math.pow(y[j]-y[j+1],2))
            if length-cellsize > slp_tol:
                x1=(length*x[j]-cellsize/2*(x[j]-x[j+1]))/length
                y1=(length*y[j]-cellsize/2*(y[j]-y[j+1]))/length
                x2=(length*x[j]-(length-cellsize/2)*(x[j]-x[j+1]))/length
                y2 = (length * y[j] - (length - cellsize / 2) * (y[j] - y[j + 1])) / length
                sx_old.append(x1)
                sx_old.append(x2)
                sy_old.append(y1)
                sy_old.append(y2)
                # sx_old=[sx_old,x1,x2]
                # sy_old=[sy_old,y1,y2]
                if j==limit:
                    x0=sx_old[newparts1[i]]
                    y0=sy_old[newparts1[i]]
                    # sx_old=[sx_old,x0]
                    # sy_old=[sy_old,y0]
                    sx_old.append(x0)
                    sy_old.append(y0)
            else:
                    if j<limit:
                        sx_old.append(float((x[j]+x[j+1])/2))
                        sy_old.append(float((y[j]+y[j+1])/2))
                        # sx_old=[sx_old,float((x[j]+x[j+1])/2)]
                        # sy_old=[sy_old,float((y[j]+y[j+1])/2)]
                    else:
                        a=sx_old[newparts1[i]]
                        b=sy_old[newparts1[i]]
                        # sx_old=[sx_old,float((x[j]+x[j+1])/2),sx_old[newparts1[i]]]
                        # sy_old=[sy_old,float((y[j]+y[j+1])/2),sy_old[newparts1[i]]]
                        sx_old.append(float((x[j]+x[j+1])/2))
                        sx_old.append(a)
                        sy_old.append(float((y[j]+y[j+1])/2))
                        sy_old.append(b)
    newparts=newparts1
    npts=len(sx_old)
    bFrist=1
    sx=[]
    sy=[]
    n=0
    for i in range(0,n_parts):

        pre_result,pre_slope=Line_Slope(sx_old[newparts1[i]],sy_old[newparts1[i]],sx_old[newparts1[i]+1],sy_old[newparts1[i]+1])

        if bFrist:
            sx.append(sx_old[newparts1[i]])
            sy.append(sy_old[newparts1[i]])
            bFrist=0
        else:
            sx.appned(sx_old[newparts1[i]])
            sy.appned(sy_old[newparts1[i]])
        if i == n_parts-1:
            limit=npts-2
        else:
            limit=newparts1[i+1]-2
        for j in range(newparts1[i]+1,limit+1):
            cur_result,cur_slope=Line_Slope(sx_old[j],sy_old[j],sx_old[j+1],sy_old[j+1])
            if pre_result !=cur_result:
                sx.append(sx_old[j])
                sy.append(sy_old[j])
                # sx=[sx,sx_old[j]]
                # sy=[sy,sy_old[j]]
                pre_result=cur_result
                pre_slope=cur_slope
            else:
                if cur_result==1:
                    pre_result=cur_result
                    pre_slope=cur_slope
                    if i<n_parts-1:
                        newparts[(i+1):(n_parts-1)]=newparts[(i+1):(n_parts-1)]-1
                    continue
                else:
                    if abs(pre_slope-cur_slope)>slp_tol:
                        sx.append(sx_old[j])
                        sy.append(sy_old[j])
                        # sx=[sx,sx_old[j]]
                        # sy=[sy,sy_old[j]]
                        pre_result=cur_result
                        pre_slope=cur_slope
                    else:
                        n = n + 1
                        if i <n_parts-1:
                            for k in range(i+1,n_parts):
                                newparts[k]=newparts[k]-1
                            # newparts[(i+1):(n_parts-1)]=newparts[(i+1):(n_parts-1)]-1
                        continue
        sx.append(sx_old[newparts1[i]])
        sy.append(sy_old[newparts1[i]])
        # sx=[sx,sx_old[newparts1[i]]]
        # sy=[sy,sy_old[newparts1[i]]]
    return 0,newparts,sx,sy
def Line_Slope(x1,y1,x2,y2):
    slope=0.0
    if x1==x2:
        return 1,slope
    slope=float(y2-y1)/float(x2-x1)
    return 0,slope

def usage():
    print("""
      usage Smoothing_shaoe.py <inputfilePath> <outputfilePath> <smoothing_order> <cell_size>

        inputfilePath  (input)  shapefile
        outputfilePath (output) shapefile
        smoothing_order (input) Number of cycles
        cell_size       (input) cell size
        Example:
       Smoothing_shape(D:\作业\冰川\处理锯齿\14_rgi60_SouthAsiaWest_test\14_rgi60_SouthAsiaWest.shp,D:\作业\冰川\处理锯齿\smoothing.shp,4,1)
    """)
    sys.exit(-1)
def main():
    inputfile=sys.argv[1]
    outputfile=sys.argv[2]
    smoothing_order=int(sys.argv[3])
    cell_size=sys.argv[4]
    if len(sys.argv)<2:
        usage()
    infile=fiona.open(inputfile,mode='r')
    schema={'geometry':'Polygon','properties':{'RGIId':'str','GLIMSId':'str','BgnDate':'str','EndDate':'str','CenLon':'float','CenLat':'float',
                                               'O1Region':'str','O2Region':'str','Area':'float','Zmin':'int','Zmax':'int','Zmed':'int',
                                               'Slope':'float','Aspect':'int','Lmax':'int','Status':'int','Connect':'int','Form':'int',
                                               'TermType':'int','Surging':'int','Linkages':'int','Name':'str'}}
    # schema=infile.schema.copy()
    outfile=fiona.open(outputfile,mode='w',driver='ESRI Shapefile',schema=schema,crs='EPSG:4326',encoding='utf-8')
    r=shapefile.Reader(inputfile)
    smoothing_shape(infile,r,smoothing_order,cell_size,outfile)



if __name__=="__main__":
    main()