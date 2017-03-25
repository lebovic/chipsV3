#!/usr/bin/env python
"""This script parses a set of _conservation.R and fills that into a
template (conserv_plotConserv.R.txt) to generate an R script 
(conserv_plotConserv.R) that, when run, generates the plot (as a png)
"""
import os
import sys
from optparse import OptionParser
from string import Template

def parseR(rscript):
    """Given a {run}_conservation.R, parse out the X,Y data--returns a tuple
    """
    f = open(rscript)
    #SKIP lines 1-3
    tmp = f.readline().strip()
    tmp = f.readline().strip()
    tmp = f.readline().strip()

    #LINE 4: READ in x -- ignored for now
    tmp = f.readline().strip()
    x = eval(tmp.replace("x<-c",""))

    #READ in y -- ignored for now
    tmp = f.readline().strip()
    y = eval(tmp.replace("y0<-c",""))

    #print(len(x), len(y))
    return (x,y)
    
def dataToString(var, data):
    """Given a tuple of data, and a name to save it as
    returns var <- c(data)
    """
    #convert data to strings
    d = [str(d) for d in data]
    return "%s <- c(%s)" % (var, ",".join(d))

def main():
    usage = "USAGE: %prog -r [R file 1] -r [R file 2] ...-r [R file N] -t [path to conserv_plotConserv.R] -o [output.png]; OUTPUTS- conserv_plotConserv.R"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-r", "--rfiles", action="append", help="list of files")
    optparser.add_option("-t", "--template", help="conserv_plotConserv.R file path")
    optparser.add_option("-p", "--pngout", help="path of png output (not produced)")
    optparser.add_option("-o", "--output", help="conserv_plotConserv.R output", default="./conserv_plotConserv.R")
    (options, args) = optparser.parse_args(sys.argv)

    if (not options.rfiles or not options.template or not options.pngout):
        optparser.print_help()
        sys.exit(-1)

    #GET run names
    runs = [r.split("/")[-1].replace("_conserv.R", "") for r in options.rfiles]
    #PARSE the data from the runs
    data = { runName: parseR(rscript) for (runName, rscript) in zip (runs,options.rfiles) }
    #compose the data string: RUN <- c(...YDATA...)
    data_str = "\n".join([dataToString(r, data[r][1]) for r in runs])
    
    #runNames = c(Run1, Run2,...)
    #runNamesQt = c("Run1", "Run2",...)
    runNames = "c(%s)" % ",".join(runs)
    runNamesQt = "c(\"%s\")" % "\",\"".join(runs)

    #for x labels, try to take the first elm's x-s
    x = [str(d) for d in data[runs[0]][0]]
    xvals = "c(%s)" % ",".join(x)

    #handle the plotting
    firstElm = runs[0]
    lines = ["lines(x,%s,type='l',col=icolor[%s])" % (r,i+2) for (i,r) in enumerate(runs[1:])]
    lines_str = "\n".join(lines)

    #PARSE the data
    #read in template
    f = open(options.template)
    tmp = Template(f.read())
    out = tmp.substitute(data=data_str, runNames=runNames, rowLabel=runNamesQt,
                         xvals=xvals, pngout=options.pngout, firstElm=firstElm,
                         rest=lines_str)
    f.close()

    #WRITE output conserv_plotConserv.R
    f = open(options.output, "w")
    f.write(out)
    f.close()

if __name__=='__main__':
    main()
