import numpy as np
import sys
# For binary file writing
import struct

def writePolyFTSDatFile(outfilename, griddim, hcell, complexdata, kspacedata, fielddata, fielddataimpart=[]):
    # Collect DIM, NPW, and nfields from other inputs
    Dim=len(griddim)
    NPW=1
    for i in range(Dim):
        NPW = NPW * griddim[i]
    nfields = len(fielddata)
    #
    f = open(outfilename,"w")
    f.write("# Format version 3\n")
    f.write("# nfields = {0}\n".format(nfields))
    f.write("# NDim = {0}\n".format(Dim))
    f.write("# PW grid = ")
    for i in range(Dim):
        f.write("{0} ".format(griddim[i]))
    f.write("\n# k-space data = {0} , complex data = {1}\n".format(int(kspacedata),int(complexdata)))

    if Dim == 1:
      if kspacedata:
        f.write("# Columns: kx indxx fielddata\n")
      else:
        f.write("# Columns: x fielddata\n")
      m=int(0)
      l=int(0)
      for ix in range(0,griddim[0]):
        l = ix
        if kspacedata and l > griddim[0]/2:
            l = l - griddim[0]
        if not kspacedata:
          xcart = hcell[0][0] * float(l) / griddim[0]
        else:
          xcart = hcell[0][0] * float(l)
        if kspacedata:
          f.write("%7.4f %d " % (xcart,l))
        else:
          f.write("%7.4f " % (xcart))

        if complexdata:
            for n in range(nfields):
                f.write("%16.10g %16.10g " % (fielddata[n][m], fielddataimpart[n][m]))
        else:
            for n in range(nfields):
                f.write("%16.10g " % (fielddata[n][m]))
        f.write("\n")
        m = m + 1

    elif Dim == 2:
      if kspacedata:
          f.write("# Columns: kx ky indxx indxy fielddata\n")
      else:
          f.write("# Columns: x y fielddata\n")
      m=int(0)
      lm=[0,0]
      for ix in range(0,griddim[0]):
        lm[0] = ix
        if kspacedata and lm[0] > griddim[0]/2:
          lm[0] = lm[0] - griddim[0]
        for iy in range(0,griddim[1]):
          lm[1] = iy
          if kspacedata and lm[1] > griddim[1]/2:
            lm[1] = lm[1] - griddim[1]

          xfrac = [float(i) for i in lm]
          if not kspacedata:
            xfrac = np.divide(xfrac,griddim)

          xcart = [0., 0.]
          for i in range(Dim):
            for j in range(Dim):
              xcart[j] = xcart[j] + hcell[i][j]*xfrac[i]
          if kspacedata:
            f.write("%7.4f %7.4f %d %d " % (xcart[0], xcart[1], lm[0], lm[1]))
          else:
            f.write("%7.4f %7.4f " % (xcart[0], xcart[1]))

          if complexdata:
              for n in range(nfields):
                  f.write("%16.10g %16.10g " % (fielddata[n][m], fielddataimpart[n][m]))
          else:
              for n in range(nfields):
                  f.write("%16.10g " % (fielddata[n][m]))
          f.write("\n")
          m = m + 1
        f.write("\n")

    elif Dim == 3:
      if kspacedata:
        f.write("# Columns: kx ky kz indxx indxy indxz fielddata\n")
      else:
        f.write("# Columns: x y z fielddata\n")
      m=int(0)
      lmn=[0,0,0]
      for ix in range(0,griddim[0]):
        lmn[0] = ix
        if kspacedata and lmn[0] > griddim[0]/2:
          lmn[0] = lmn[0] - griddim[0]
        for iy in range(0,griddim[1]):
          lmn[1] = iy
          if kspacedata and lmn[1] > griddim[1]/2:
            lmn[1] = lmn[1] - griddim[1]
          for iz in range(0,griddim[2]):
            lmn[2] = iz
            if kspacedata and lmn[2] > griddim[2]/2:
              lmn[2] = lmn[2] - griddim[2]

            xfrac = [float(i) for i in lmn]
            if not kspacedata:
              xfrac = np.divide(xfrac,griddim)

            xcart = [0., 0., 0.]
            for i in range(Dim):
              for j in range(Dim):
                xcart[j] = xcart[j] + hcell[i][j]*xfrac[i]
            if kspacedata:
              f.write("%7.4f %7.4f %7.4f %d %d %d " % (xcart[0], xcart[1], xcart[2], lmn[0], lmn[1], lmn[2]))
            else:
              f.write("%7.4f %7.4f %7.4f " % (xcart[0], xcart[1], xcart[2]))

            if complexdata:
                for n in range(nfields):
                    f.write("%16.10g %16.10g " % (fielddata[n][m], fielddataimpart[n][m]))
            else:
                for n in range(nfields):
                    f.write("%16.10g " % (fielddata[n][m]))
            f.write("\n")
            m = m + 1
        f.write("\n")
      f.close()

def writePolyFTSBinFile(outfilename, griddim, hcell, complexdata, kspacedata, fielddata, fielddataimpart=[]):
    # Collect DIM, NPW, and nfields from other inputs
    Dim=len(griddim)
    NPW=1
    for i in range(Dim):
        NPW = NPW * griddim[i]
    nfields = len(fielddata)
    #
    f = open(outfilename,"wb")
    f.write(struct.pack('@8s',"FieldBin"))
    f.write(struct.pack('@b',0)) # Null termination character for string
    f.write(struct.pack('@I',51)) # Version code
    f.write(struct.pack('@I',nfields))
    f.write(struct.pack('@I',Dim))
    f.write(struct.pack('@{0}L'.format(Dim),*griddim))
    f.write(struct.pack('@2?',kspacedata,complexdata))
    harray=np.reshape(hcell,(Dim*Dim))
    f.write(struct.pack('@{0}d'.format(Dim*Dim),*harray))
    # Determine element size and write - always DP here
    if complexdata:
        f.write(struct.pack('@L',16))
    else:
        f.write(struct.pack('@L',8))
    # Write field data
    if complexdata:
        f.write(struct.pack("@{0}d".format(NPW*2*nfields),*(np.reshape(np.stack((fielddata,fielddataimpart),axis=2),(NPW*2*nfields)))))
    else:
        f.write(struct.pack("@{0}d".format(NPW*nfields),*(np.reshape(fielddata,(NPW*nfields)))))
    f.close()
