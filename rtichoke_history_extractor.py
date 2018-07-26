from argparse import ArgumentParser
import io,re

parser = ArgumentParser(description="rtichoke history extractor")
parser.add_argument("-i","--inputFile",default=".rtichoke_history")
parser.add_argument("-o","--outputFile",default=".Rhistory_rtichoke")

args = parser.parse_args()
with io.open(args.inputFile,"r") as inputfile:
    with io.open(args.outputFile,"w") as outputfile:
        rhis = inputfile.readlines()
        Rhis = []
        for line in rhis:
            if line=='':
                continue
            elif line[0]=='#':
                continue
            elif line[0]=='+':
                Rhis.append(line[1:])
        outputfile.writelines(Rhis)
print("Finished!")





