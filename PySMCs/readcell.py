"""
Read SMC grid cell arrays from text files.

Function readcell(filestring) read the cell file specified by
the filestring, which includes path/file_name.

If Arctic cell file is also needed, the filestring will be a
list of both global and Arctic cell files.  It will merge the 
two parts together.

Function readtext(filestring) can also be used to read single 
cell file though it is designed to read any formated text file.
It allows descriptive lines at the beginning of file to be
skipped and uses the last skipped line as header parameters. 
Also note readtext return split header parameters as a list
while readcell returns unsplit header line or two header lines
if Arctic cell is also read.

The main() function demonstrates how the readcell could be called
as a stand alone function and merge the Arctic part manually.

              JGLi18Feb2019 

"""

##  Assume to skip the first line [0]. If more than one
##  lines to skip, specifiy the line indexes in the
##  parameter list skiprows=[0,1], for first two lines.
##  But only the first line is returned as a list.

def readtext(textfile, skiprows=[0]):
    import numpy  as np
    import pandas as pd

    print( " Read table from ", textfile)
    archd=open(textfile, 'r') 
    for i in range(len(skiprows)):
        hdrskp=archd.readline()

##  Split the last one of skipped lines as header parameters.
    hdlist = hdrskp.split()
    archd.close()

    arcel=pd.read_csv(textfile, sep='\s+',skiprows=skiprows,header=None)
    table=arcel.values

    return hdlist, table

##  Input celfiles as a list even if there is only one file.
##  For instance celfiles=['path/celfile.dat', 'path/arcfile.dat']
def readcell(celfiles):
    import numpy  as np
    import pandas as pd
    
    nfls = 0 
    for celfile in celfiles:
        print( " Read cel from ", celfile)
        archd=open(celfile, 'r') 
        hdlin=archd.readline()
        archd.close()
#   hdlist=archd.readline().split()
#   ncs=np.array( hdlist ).astype(np.int64)
         
        flcel=pd.read_csv(celfile, sep='\s+',skiprows=[0],header=None)
        celin=flcel.values
        nfls += 1

        if( nfls <= 1):
            headrs = [ hdlin ]
            cel = celin.copy()
        else:
            headrs.append( hdlin )
            cel = np.vstack( (cel, celin) )
        
    return headrs, cel


## 
def main():

    import numpy  as np
    import pandas as pd

    cellfile = input(" Please enter the cell file as 'path/file.dat' > ")
#   cellfile='./G50SMCels.dat'
    ncs, cel = readcell( [cellfile] )
    nc = ncs[0]

    Arctic = input(" Is there an Arctic part cell file? 'Yes' or 'No' > ")
#   Arctic = 'Yes'
    if Arctic.strip()[0] == 'Y':
        Arcfile = input(" Please enter the Arctic part cell file > ")
#       Arcfile='./G50SMCBAr.dat'
        nab, acel = readcell( [Arcfile] )
        nc = ncs[0] + nab[0]
        cel = np.vstack( (cel, acel) )
        print( "Arctic cells appended to global ones.") 

    print( "Total number of cells nc =", nc)
    print( "First cell array is ", cel[0,:])
    print( " Last cell array is ", cel[nc-1,:])


if __name__ == '__main__':
    main()

