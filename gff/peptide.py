'''
Descripttion: 根据蛋白质结构域的注释，获取其在基因组中的真实坐标
version: 
Author: zpliu
Date: 2022-07-25 11:16:42
LastEditors: zpliu
LastEditTime: 2022-07-25 11:18:47
@param: 
'''
#----------------------------------------
#todo 使用数据库文件
#----------------------------------------
import gffutils
import pandas as pd 
db = gffutils.FeatureDB(
    '/data/cotton/zhenpingliu/genome_data/Ghirsutum_genome_HAU_v1.1/Ghirsutum_gene_model.db',
     keep_order=True)

queryItem=('Ghir_A01G015600','cds.Ghir_A01G015600.1',40,60)
queryList=[queryItem] 

out=[]
for queryItem in queryList:
    gene=db[queryItem[0]]
    cdsStart=(queryItem[2]-1)*3+1
    cdsEnd=(queryItem[3])*3
    CDSarray=[]
    for i in db.children(gene, featuretype='CDS', order_by='start'):
        if i.attributes['ID'][0]==queryItem[1]:
            CDSarray.append(
                (i.start,i.end,abs(i.start-i.end)+1,i.strand,i.attributes['ID'][0])
            ) 
        else:
            pass 
    #* 判断正负链,
    totalSeq=0
    if CDSarray[0][3]=="+":
        #TODO 将CDS片段进行累加，判断domain的起始位点是否位于该exon上
        for cds_seq in CDSarray:
            totalSeq+= cds_seq[2]
            if cdsStart>totalSeq:
                #* 还没有走到对应的exon
                pass 
            else:
                realStart= cds_seq[1]-(totalSeq-cdsStart)
                totalSeq=0
                break 
        
        for cds_seq in CDSarray:
            totalSeq+= cds_seq[2]
            if totalSeq >=cdsEnd:
                #TODO: 走到了对应的exon
                realEnd=cds_seq[1]-(totalSeq-cdsEnd)
                totalSeq=0
                break
            else:
                pass 
    else:
        for i in range(len(CDSarray)-1,-1,-1):
            totalSeq+= CDSarray[i][2]
            if cdsStart>totalSeq:
                pass 
            else:
                realStart=CDSarray[i][0]+(totalSeq-cdsStart)
                totalSeq=0
                break 
        for i in range(len(CDSarray)-1,-1,-1):
            totalSeq+= cds_seq[2]
            if totalSeq >=cdsEnd:
                realEnd= CDSarray[i][0]+(totalSeq-cdsEnd)
                totalSeq=0
                break
    # #------------------------------------
    # #* 将真实坐标的起始与所有的cds片段进行交集
    # #------------------------------------
    intersectArray=[]
    for cds_seq in CDSarray:
        if  cds_seq[0]>min(realStart,realEnd) and  cds_seq[0]<max(realStart,realEnd):
            intersectArray.append(cds_seq[0])
        if  cds_seq[1]>min(realStart,realEnd) and  cds_seq[1]<max(realStart,realEnd):
            intersectArray.append(cds_seq[1])
    intersectArray.append(realStart)
    intersectArray.append(realEnd)
    intersectArray.sort()
    intersectExon=[]
    for i in range(0,len(intersectArray),2):
        intersectExon.append(
            "{}-{}".format(intersectArray[i],intersectArray[i+1])
        )
    out.append(
        (queryItem[0],queryItem[1],queryItem[2],queryItem[3],gene.seqid,";".join(intersectExon),gene.strand)
    )
    
out=pd.DataFrame(out,columns=[
    'geneId','cdsId','peptideStart','peptideEnd','chrom','genomeRegion','strand'
])    