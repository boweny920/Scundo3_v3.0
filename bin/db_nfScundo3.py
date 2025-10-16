#!/usr/bin/env python

import os
import json
import csv
import sys


with open(sys.argv[1], 'r') as f:
    data = f.readlines() 
    f.close()

orders = []
empty = " "
noorderType=0
noreadType=0
noreadLength=0
nogenome=0
nogenomes=[]
noanalysisGoals=0

for index in data:
    order = {}
    i = index.strip('\n').rstrip("/")
    order['Order'] = i.split('/')[-1]
    order['Lab'] = i.split('/')[-3]
    order['Requester'] = i.split('/')[-2]

    if os.path.exists(i+'/targets.tsv') != True:
        continue
    
    with open(i+'/targets.tsv') as tsvfile:
        reader = csv.DictReader(tsvfile, dialect='excel-tab')
        
        '''For in-houst secundo'''
        if "Flowcell" in reader.fieldnames:
            reader = sorted(reader, key=lambda d:d['Secundo Name'])
            record = []
            simples = []
            simple = {}
            simple['Flowcell'] = []
            column = [row for row in reader]
            countfor = 0
            # initid = column[0]['Library']
            initid = ""
            initflowid = ""
            for row in column:
                first = row['Secundo Name']
                if(first == initid ): 
                    # go add the flowcell
                    firstflowid = row['Flowcell']
                    if (initflowid == firstflowid):
                        continue
                    else:
                        simple['Flowcell'].append(row['Flowcell'])
                        initflowid = row['Flowcell']
                else:
                    # create a new simple.
                    simple = {}
                    simple['Flowcell']=[]   
                    simple['LibraryID'] = row['Library']
                    simple['SecundoName'] = row['Secundo Name']
                    simple['SampleName'] = row['Sample Name']
                    simple['Flowcell'].append(row['Flowcell'])
                    initflowid = row['Flowcell']
                    initid = first
                    simples.append(simple)
            order['Samples']=simples
            tsvfile.close()
        else:
            """For pubdata roboSecundo"""
            reader = sorted(reader, key=lambda d:d['Secundo Name'])
            record = []
            simples = []
            simple = {}
            simple['Flowcell'] = []
            column = [row for row in reader]
            countfor = 0
            # initid = column[0]['Library']
            initid = ""
            initflowid = ""
            for row in column:
                first = row['Secundo Name']
                if(first == initid ): 
                    # go add the flowcell
                    firstflowid = row['Order']
                    if (initflowid == firstflowid):
                        continue
                    else:
                        simple['Flowcell'].append(row['Order'])
                        initflowid = row['Order']
                else:
                    # create a new simple.
                    simple = {}
                    simple['Flowcell']=[]   
                    simple['LibraryID'] = row['Library']
                    simple['SecundoName'] = row['Secundo Name']
                    simple['SampleName'] = row['Sample Name']
                    simple['Flowcell'].append(row['Order'])
                    initflowid = row['Order']
                    initid = first
                    simples.append(simple)
            order['Samples']=simples
            tsvfile.close()
    
    with open(i+'/lims_order.json','r') as jsonfile:
        load_dict = json.load(jsonfile) 
        if isinstance(load_dict,list):
            load_dict = load_dict[0]
        else: 
            load_dict = load_dict
            

        # fetch data
        if ('orderType' in load_dict):
            order['OrderType']=load_dict['orderType']
        else:
            order['OrderType']=empty
            noorderType += 1
        
        if ('readType' in load_dict):
            order['ReadType']=load_dict['readType']
        else:
            order['ReadType']=empty
            noreadType += 1
        
        if ('readLength' in load_dict):
            order['ReadLength']=load_dict['readLength']
        else:
            order['OrderType']=empty
            noreadLength += 1
        
        if ('genome' in load_dict):
            order['Genome']=load_dict['genome']
        else:
            order['Genome']=empty
            nogenome += 1
            nogenomes.append(order['Order'])
        
        if ('analysisGoals' in load_dict):
            order['AnalysisGoals']=load_dict['analysisGoals']
        else:
            order['AnalysisGoals']=empty
            noanalysisGoals += 1
            

        jsonfile.close()  

    orders.append(order)

j = json.dumps(orders)

with open("/n/core/Bioinformatics/tools/biotools/db/datafile.json","w") as dump_f:
    json.dump(orders,dump_f)