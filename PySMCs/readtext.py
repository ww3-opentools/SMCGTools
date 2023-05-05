"""
Read data arrays from a tabled text file. Missing values in the
table will be filled as none.

Function readtext(filestring) read the cell file specified by
the filestring, which includes path/file_name.txt

The main() function demonstration how the readtext could be called
as a stand alone function.

              JGLi18Feb2019 

"""

##  Assume to skip the first line [0]. If more than one
##  lines to skip, specifiy the line indexes in the
##  parameter list skiprows=[0,1], for first two lines.
##  But only the last skipped line is returned as a list.

def readtext(textfile, skiprows=[0]):
    import numpy  as np
    import pandas as pd

    print (" Read table from ", textfile)
    archd=open(textfile, 'r') 
##  Assume last line in skiprows contains the header parameters.
    for i in range(len(skiprows)):
        hdrskp=archd.readline()

##  Split the last one of skipped lines as header parameters.
    hdlist = hdrskp.split()
    archd.close()

    arcel=pd.read_csv(textfile, sep='\s+',skiprows=skiprows,header=None)
    table=arcel.values

    return hdlist, table

def main():

    import numpy  as np
    import pandas as pd

    cellfile = input(" Please enter the cell file as 'path/file.dat' > ")
#   cellfile='./G50SMCels.dat'
    ncs, cel = readtext( cellfile )
    nc = int(ncs[0])

    print ("Total number of cells nc =", nc)
    print ("First cell array is ", cel[0,:])
    print (" Last cell array is ", cel[nc-1,:])

if __name__ == '__main__':
    main()

