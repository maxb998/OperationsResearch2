import os, argparse, re, csv, glob
import numpy as np


class InputParams:
    inputDir = ''
    prefix = ''
    outFname = ''
    separator = ''
    lineList = []
    stopLineList = []

    def __init__(self, inputDir:str, prefix:str, outFname:str, separator:str, lineList:list, stopLineList:list):
        self.inputDir = inputDir
        self.prefix = prefix
        self.outFname = outFname
        self.separator = separator
        self.lineList = lineList
        self.stopLineList = stopLineList

class FilenameExtractedData:
    tuningVar = ''
    instanceName = ''
    timeLimit = ''

    def __init__(self, tuningVar:str, instanceName:str, timeLimit:str):
        self.tuningVar = tuningVar
        self.instanceName = instanceName
        self.timeLimit = timeLimit


def argParser() -> InputParams:
    parser = argparse.ArgumentParser(prog='Log2csv', epilog='Log 2 Csv', description='Read data from multiple log files and outputs a csv file ready for performance profile')
    parser.add_argument('-I', '--inputDir', metavar='str', required=True, type=str, help='Directory containing all the log files to read')
    parser.add_argument('-P', '--prefix', metavar='str', required=True, type=str, help='Prefix present in all the files inside the directory')
    parser.add_argument('-O', '--outputFname', metavar='str', required=True, type=str, help='Output filename for the IterCount csv')
    parser.add_argument('-S', '--separator', choices=[';','space',','], default=';', required=False, type=str, help='Type of separator for csv file')
    parser.add_argument('-L', '--lineNames', required=True, nargs='+', default=['Final cost = '], type=str, help='String contained in the line/lines that contain the number to extract')
    parser.add_argument('--stopLines', required=False, nargs='+', default=[], type=str, help='String contained in the line/lines that will stop the reading of the file')

    args = parser.parse_args()

    if not os.path.isdir(args.inputDir):
        print('Directory specified as --inputDir is not a directory')
        exit()
    
    if not args.outputFname.endswith('.csv'):
        args.outputFname += '.csv'

    if os.path.isfile(args.outputFname):
        print('Output file already exists')
        exit()
    
    if args.separator == 'space':
        args.separator = ' '

    return InputParams(args.inputDir, args.prefix, args.outputFname, args.separator, args.lineNames, args.stopLines)

def getInfoFromFilename(filename:str, params:InputParams) -> FilenameExtractedData:
    basename = os.path.basename(filename)

    if params.prefix not in basename:
        print('File \"' + basename + '\" does not contain the prefix \"' + params.prefix + '\"')
        exit()
    
    fnameSplitted = basename.split('_')
    compsInPrefix = len(params.prefix.split('_'))

    return FilenameExtractedData(
        instanceName=fnameSplitted[compsInPrefix + 0],
        timeLimit=fnameSplitted[compsInPrefix + 1][:len(fnameSplitted[compsInPrefix + 1]) - 1],
        tuningVar=fnameSplitted[compsInPrefix + 2]
    )
    
def getFilesWithSameInstanceAndTimeLimit(filelist:list) -> [str, ...]:
    basename = os.path.basename(filelist[0])
    prefix = basename[0:basename.rfind('_')]
    prefix = prefix[0:prefix.rfind('_')] + '_'

    retVal = []
    for i in range(len(filelist)):
        if prefix in os.path.basename(filelist[i - len(retVal)]):
            retVal.append(filelist.pop(i - len(retVal)))

    return retVal

def getFilesWithDifferentSeedRuns(flist:list) -> [str, ...]:
    basename = os.path.basename(flist[0])
    prefix = basename[0:basename.rfind('_')] + '_'

    retVal = []
    for i in range(len(flist)):
        if prefix in os.path.basename(flist[i - len(retVal)]):
            retVal.append(flist.pop(i - len(retVal)))

    return retVal

def readLogFile(filename:str, params:InputParams) -> float:
    logFile = open(filename)

    stopRead = False
    retVal = -1.0
    for line in logFile:
        for lineName in params.lineList:
            if lineName in line:
                retVal = float(re.findall(r'[-+]?(?:\d*\.*\d+)', line[line.find(lineName) + len(lineName):])[0])
        for stopLine in params.stopLineList:
            if stopLine in line:
                stopRead = True
        if stopRead:
            break

    logFile.close()
    
    if retVal == -1.0:
        print('Could not find keystring in \"' + filename + '\"')
        exit()
    
    return retVal

def main():
    params = argParser()
    
    filelist = sorted(glob.glob(os.path.join(params.inputDir, '*.log')))

    for f in filelist:
        if params.prefix not in f:
            print('File \"' + f + '\" is missing prefix \"' + params.prefix + '\"')

    # Sort by number of nodes
    nNodes = np.empty(shape=len(filelist), dtype=np.int32)
    for i in range(len(filelist)):
        nNodes[i] = re.findall(r'\d+', filelist[i])[0]
    sortedArgs = nNodes.argsort()
    filelist = [filelist[i] for i in sortedArgs]
    nNodes = [nNodes[i] for i in sortedArgs]
    del sortedArgs

    # Subsort for each number of nodes by the name of the instance
    i = 0
    while i < len(nNodes):
        j = i
        while i < len(nNodes) and  nNodes[j] == nNodes[i]:
            i += 1
        filelist[j:i] = sorted(filelist[j:i], key=lambda p: getInfoFromFilename(p,params).instanceName )
    del nNodes
    
    # Subsort for each instance by time limit
    i = 0
    while i < len(filelist):
        j = i
        while i < len(filelist) and getInfoFromFilename(filelist[i], params).instanceName == getInfoFromFilename(filelist[j], params).instanceName:
            i += 1
        filelist[j:i] = sorted(filelist[j:i], key=lambda p: float(getInfoFromFilename(p,params).timeLimit))

    # Count number of columns to write into the csv
    tuningVarsList = []
    for f in filelist:
        fnameParams = getInfoFromFilename(os.path.basename(f), params)
        if fnameParams.tuningVar not in tuningVarsList:
            tuningVarsList.append(fnameParams.tuningVar)

    tuningVarsList.sort(key=lambda i: -float(i))

    print('Detected tuning vars:')
    print(tuningVarsList)

    # for f in filelist:
    #     print('  ' + os.path.basename(f))

    csvfile = open(params.outFname, 'a')

    csvwriter = csv.writer(csvfile, delimiter=params.separator)

    csvline = [str(len(tuningVarsList))]
    for var in tuningVarsList:
        csvline.append(var)
    csvline.append('TimeLimit')
    csvwriter.writerow(csvline)

    while len(filelist) > 0:

        filesWithSameInstAndTimeLimit = getFilesWithSameInstanceAndTimeLimit(filelist)
        # print('  ' + '\n  '.join(filesWithSameInstAndTimeLimit))
        # print()

        fnameInfo = getInfoFromFilename(filesWithSameInstAndTimeLimit[0], params)

        csvline[0] = fnameInfo.instanceName
        for i in range(1, len(csvline)-1, 1):
            csvline[i] = 'null'
        csvline[len(csvline)-1] = fnameInfo.timeLimit

        while len(filesWithSameInstAndTimeLimit) > 0:
            
            dataList = []

            filesWithDifferentSeed = getFilesWithDifferentSeedRuns(filesWithSameInstAndTimeLimit)
            #print('\t' + '\n\t'.join(filesWithDifferentSeed))

            for i in range(len(filesWithDifferentSeed)):
                dataList.append(readLogFile(filesWithDifferentSeed[i], params))

            #print(dataList)

            avgStr = str(sum(dataList)/len(dataList))

            fnameInfo = getInfoFromFilename(filesWithDifferentSeed[0], params)

            for i in range(len(tuningVarsList)):
                if fnameInfo.tuningVar == tuningVarsList[i]:
                    csvline[i+1] = avgStr

        print(csvline)
        csvwriter.writerow(csvline)
    
    csvfile.close()
   

if __name__ == '__main__' :
    main()
