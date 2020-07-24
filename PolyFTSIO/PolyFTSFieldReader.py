import numpy as np
import sys
# For binary file reading
import struct
# gzip support
import gzip

class PolyFTSFieldReader(object):
    """ A class to read, store and manupulate data from PolyFTS field files """

    def __init__(self):
        self._gzipped = False
        self._BinFile = False


    def _testBinFile(self,filepath):
        self._BinFile = False
        if self._gzipped:
            i = gzip.open(filepath,"rb")
        else:
            i = open(filepath,"rb")
        header = struct.unpack_from("@8s",i.read(8))[0] # Check first 8 bytes
        self.version = struct.unpack_from("@I",i.read(5),offset=1)[0] # Note the need to skip one byte due to the null termination character that got saved in the C string
        i.close()
        if header == "FieldBin" and self.version == 51:
            self._BinFile = True


    def _testgzipFile(self,filepath):
        self._gzipped = False
        ftest = open(filepath,'rb')
        twobytes=struct.unpack('2c',ftest.read(2)) # read 2 bytes and unpack as a tuple of chars
        ftest.close()
        if twobytes[0] == '\x1f' and twobytes[1] == '\x8b':
            self._gzipped = True


    def _generateMeshCoordinates_rspace(self):
        self.AllMeshCoords = np.zeros((self.Dim,self.NPW))

        if self.Dim  == 1:
            m=int(0)
            l=int(0)
            for ix in range(0,self.griddim[0]):
                l = ix
                xcart = self.hcell[0][0] * float(l) / self.griddim[0]

                self.AllMeshCoords[:,m] = xcart
                m+=1

        elif self.Dim == 2:
            m=int(0)
            lm=[0,0]
            for ix in range(0,self.griddim[0]):
                lm[0] = ix
                for iy in range(0,self.griddim[1]):
                    lm[1] = iy

                    xfrac = [float(i) for i in lm]
                    xfrac = np.divide(xfrac,self.griddim)

                    xcart = [0., 0.]
                    for i in range(self.Dim):
                        for j in range(self.Dim):
                            xcart[j] = xcart[j] + self.hcell[i][j]*xfrac[i]
                    self.AllMeshCoords[:,m] = xcart
                    m+=1

        elif self.Dim == 3:
            m=int(0)
            lmn=[0,0,0]
            for ix in range(0,self.griddim[0]):
                lmn[0] = ix
                for iy in range(0,self.griddim[1]):
                    lmn[1] = iy
                    for iz in range(0,self.griddim[2]):
                        lmn[2] = iz
                        xfrac = [float(i) for i in lmn]
                        xfrac = np.divide(xfrac,self.griddim)

                        xcart = [0., 0., 0.]
                        for i in range(self.Dim):
                            for j in range(self.Dim):
                                xcart[j] = xcart[j] + self.hcell[i][j]*xfrac[i]
                        self.AllMeshCoords[:,m] = xcart
                        m+=1

    def _generateMeshCoordinates_kspace(self):
        self.AllMeshCoords = np.zeros((self.Dim,self.NPW))

        if self.Dim  == 1:
            m=int(0)
            l=int(0)
            for ix in range(0,self.griddim[0]):
                l = ix
                if l > self.griddim[0]/2:
                    l = l - self.griddim[0]

                xcart = self.hcell[0][0] * l
                self.AllMeshCoords[:,m] = xcart
                m+=1

        elif self.Dim == 2:
            m=int(0)
            lm=[0,0]
            for ix in range(0,self.griddim[0]):
                lm[0] = ix
                if lm[0] > self.griddim[0]/2:
                    lm[0] = lm[0] - self.griddim[0]
                for iy in range(0,self.griddim[1]):
                    lm[1] = iy
                    if lm[1] > self.griddim[1]/2:
                        lm[1] = lm[1] - self.griddim[1]

                    xcart = [0., 0.]
                    for i in range(self.Dim):
                        for j in range(self.Dim):
                            xcart[j] = xcart[j] + self.hcell[i][j]*lm[i]
                    self.AllMeshCoords[:,m] = xcart
                    m+=1

        elif self.Dim == 3:
            m=int(0)
            lmn=[0,0,0]
            for ix in range(0,self.griddim[0]):
                lmn[0] = ix
                if lmn[0] > self.griddim[0]/2:
                    lmn[0] = lmn[0] - self.griddim[0]
                for iy in range(0,self.griddim[1]):
                    lmn[1] = iy
                    if lmn[1] > self.griddim[1]/2:
                        lmn[1] = lmn[1] - self.griddim[1]
                    for iz in range(0,self.griddim[2]):
                        lmn[2] = iz
                        if lmn[2] > self.griddim[2]/2:
                            lmn[2] = lmn[2] - self.griddim[2]

                        xcart = [0., 0., 0.]
                        for i in range(self.Dim):
                            for j in range(self.Dim):
                                xcart[j] = xcart[j] + self.hcell[i][j]*lmn[i]
                        self.AllMeshCoords[:,m] = xcart
                        m+=1



    def _readBinFile(self, filepath):
        """ Wrapper function to parse field grid data from PolyFTS unformatted binary file"""
        if self._gzipped:
            i = gzip.open(filepath,'rb')
        else:
            i = open(filepath,'rb')
        # Read the entire file
        contents = i.read()
        # Close file
        i.close()

        # Check header and file version number
        header = struct.unpack_from("@8s",contents)
        if header[0] != "FieldBin":
            sys.stderr.write("\nError: Not an unformatted Field file\n")
            sys.exit(1)
        pos = 9
        self.version = struct.unpack_from("@I",contents,offset=pos)
        pos = pos + 4
        if self.version[0] != 51:
            sys.stderr.write("\nError: Only version code 51 is currently supported\n")
            sys.exit(1)

        self.nfields = struct.unpack_from("@I",contents,offset=pos)[0]
        pos = pos + 4

        self.Dim = struct.unpack_from("@I",contents,offset=pos)[0]
        pos = pos + 4

        self.griddim = struct.unpack_from("@{0}L".format(self.Dim),contents,offset=pos)
        pos = pos + 8*self.Dim
        self.NPW = 1
        for i in range(self.Dim):
            self.NPW = self.NPW * self.griddim[i]

        (self.kspacedata,self.complexdata) = struct.unpack_from("@2?",contents,offset=pos)
        pos = pos + 2

        harray = struct.unpack_from("@{0}d".format(self.Dim*self.Dim),contents,offset=pos)
        pos = pos + 8*self.Dim*self.Dim
        self.hcell = np.reshape(harray,(self.Dim,self.Dim))

        # Determine the number of bytes per element
        elemsize = struct.unpack_from("@L",contents,offset=pos)[0]
        pos = pos + 8
        if elemsize == 4 and not self.complexdata:
            fielddata = struct.unpack_from("@{0}f".format(self.NPW*self.nfields),contents,offset=pos)
        elif elemsize == 8 and self.complexdata:
            fielddata = struct.unpack_from("@{0}f".format(2*self.NPW*self.nfields),contents,offset=pos)
        elif elemsize == 8 and not self.complexdata:
            fielddata = struct.unpack_from("@{0}d".format(self.NPW*self.nfields),contents,offset=pos)
        elif elemsize == 16 and self.complexdata:
            fielddata = struct.unpack_from("@{0}d".format(2*self.NPW*self.nfields),contents,offset=pos)
        else:
            sys.stderr.write("\nError: Unknown element size")
            sys.exit(1)

        # Now we need to return a numpy array with two indices:
        #  1 = fieldidx (with imaginary part as a distinct field index)
        #  2 = PW indx
        # For complex fields, currently the real/imaginary part is the fastest index,
        # so we will have to do a selective transpose
        # before reshaping the re/im sequence into the field indices
        if self.complexdata:
            # Copy real and imaginary parts to separate arrays
            self.AllFields       = np.array(fielddata).reshape([self.nfields,self.NPW,2]).transpose()[0].transpose()
            self.AllFieldsImPart = np.array(fielddata).reshape([self.nfields,self.NPW,2]).transpose()[1].transpose()
            # OLD way to interleave Re and Im fields
            #self.AllFields = np.array(fielddata).reshape([self.nfields,self.NPW,2]).transpose((0,2,1)).reshape([self.nfields*2,self.NPW])
        else:
            self.AllFields = np.array(fielddata).reshape([self.nfields,self.NPW])

        # Since the bin file doesn't contain the mesh coordinates, we need to generate them from the box tensor and griddim
        if self.kspacedata:
            self._generateMeshCoordinates_kspace()
        else:
            self._generateMeshCoordinates_rspace()
        # Test whether the cell is orthorhombic
        self.orthorhombic = True
        for i in range(self.Dim):
            for j in range(self.Dim):
                if i==j:
                    continue
                if abs(self.hcell[i][j]) > 1e-8:
                    self.orthorhombic = False
        # Store grid spacing
        self.spacing = [None]*self.Dim
        if self.Dim == 1:
            self.spacing[0] = self.AllMeshCoords[0][1]
        if self.Dim == 2:
            self.spacing[0] = self.AllMeshCoords[0][self.griddim[1]]
            self.spacing[1] = self.AllMeshCoords[1][1]
        if self.Dim == 3:
            self.spacing[0] = self.AllMeshCoords[0][self.griddim[1]*self.griddim[2]]
            self.spacing[1] = self.AllMeshCoords[1][self.griddim[2]]
            self.spacing[2] = self.AllMeshCoords[2][1]
        # DONE


    def _readDatFile(self, filepath):
        """ Wrapper function to parse field grid data from PolyFTS formatted file"""
        if self._gzipped:
            f = gzip.open(filepath,'r')
        else:
            f = open(filepath)
        # Now parse the header to obtain the grid dimension and other format information
        self.version=int(f.readline().strip().split()[3])
        self.nfields=int(f.readline().strip().split()[3])
        self.Dim=int(f.readline().strip().split()[3])
        self.griddim=f.readline().strip().split()[4:4+self.Dim]
        self.griddim=[int(i) for i in self.griddim] # Convert all list entries to int
        self.NPW = 1
        for i in range(self.Dim):
            self.NPW = self.NPW * self.griddim[i]
        self.kspacedata=True
        self.complexdata=True
        flagsline=f.readline().strip().split()
        if flagsline[4] == "0":
            self.kspacedata=False
        if flagsline[9] == "0":
            self.complexdata=False
        f.readline() # Skip the last header line

        # Get the field data and grid coordinates
        AllData = np.loadtxt(f).transpose()
        self.AllMeshCoords = AllData[:self.Dim]
        if not self.kspacedata:
            TMP = AllData[self.Dim:]
        else:
            TMP = AllData[2*self.Dim:] # Extra coordinate columns for kspace data
        if self.complexdata:
            self.AllFields       = TMP.reshape([self.nfields,2,self.NPW])[:,0,:]
            self.AllFieldsImPart = TMP.reshape([self.nfields,2,self.NPW])[:,1,:]
        else:
            self.AllFields       = TMP.reshape([self.nfields,self.NPW])[:,:]
        # Double the field count if we have imaginary parts

        # TODO: future versions of dat file should store the full box tensor,
        # Here we recompute it from the mesh, but there can be a significant loss of accuracy due to
        # truncation of the precision in the .dat file mesh.
        # Check whether the cell is orthorhombic and store grid spacing
        self.orthorhombic = True
        self.spacing = [None]*self.Dim
        self.hcell = np.reshape([None]*self.Dim*self.Dim,(self.Dim,self.Dim))
        if self.Dim == 1:
            self.spacing[0] = self.AllMeshCoords[0][1]
            # Reconstruct cell tensor
            if self.kspacedata == True:
                self.hcell[0][0] = self.AllMeshCoords[0][1] * self.griddim[0]
            else:
                self.hcell[0][0] = self.AllMeshCoords[0][1]
        if self.Dim == 2:
            if self.AllMeshCoords[0][1] != 0.:
                self.orthorhombic = False
            if self.AllMeshCoords[1][self.griddim[1]] != 0.:
                self.orthorhombic = False
            self.spacing[0] = self.AllMeshCoords[0][self.griddim[1]]
            self.spacing[1] = self.AllMeshCoords[1][1]
            # Reconstruct cell tensor
            if self.kspacedata == True:
                for i in range(self.Dim):
                    self.hcell[0][i] = self.AllMeshCoords[i][self.griddim[1]]
                    self.hcell[1][i] = self.AllMeshCoords[i][1]
            else:
                for i in range(self.Dim):
                    self.hcell[0][i] = self.AllMeshCoords[i][self.griddim[1]] * self.griddim[0]
                    self.hcell[1][i] = self.AllMeshCoords[i][1] * self.griddim[1]
        if self.Dim == 3:
            # Check that the first non-zero mesh point has only a z translation
            if self.AllMeshCoords[0][1] != 0. or self.AllMeshCoords[1][1] != 0.:
                self.orthorhombic = False
            # Check that the first y translation has no x,z
            if self.AllMeshCoords[0][self.griddim[2]] != 0 or self.AllMeshCoords[2][self.griddim[2]] != 0:
                self.orthorhombic = False
            # Check that the first x translation has no y,z
            if self.AllMeshCoords[1][self.griddim[2]*self.griddim[1]] != 0 or self.AllMeshCoords[2][self.griddim[2]*self.griddim[1]] != 0:
                self.orthorhombic = False
            self.spacing[0] = self.AllMeshCoords[0][self.griddim[1]*self.griddim[2]]
            self.spacing[1] = self.AllMeshCoords[1][self.griddim[2]]
            self.spacing[2] = self.AllMeshCoords[2][1]
            # Reconstruct cell tensor
            if self.kspacedata == True:
                for i in range(self.Dim):
                    self.hcell[0][i] = self.AllMeshCoords[i][self.griddim[1]*self.griddim[2]]
                    self.hcell[1][i] = self.AllMeshCoords[i][self.griddim[2]]
                    self.hcell[2][i] = self.AllMeshCoords[i][1]
            else:
                for i in range(self.Dim):
                    self.hcell[0][i] = self.AllMeshCoords[i][self.griddim[1]*self.griddim[2]] * self.griddim[0]
                    self.hcell[1][i] = self.AllMeshCoords[i][self.griddim[2]] * self.griddim[1]
                    self.hcell[2][i] = self.AllMeshCoords[i][1] * self.griddim[2]


        ## Get the grid data
        #allrows = np.loadtxt(f,usecols=range(0,ndim)+[col])
        #allcols = np.transpose(allrows)
        #gridcoords = allcols[0:ndim]
        #Lmin = [1.0 * np.min(L) / N * (N+1) for L, N in zip(gridcoords, griddim)]
        #Lmax = [1.0 * np.max(L) / N * (N+1) for L, N in zip(gridcoords, griddim)]
        #tmp = np.reshape(zip(Lmin,Lmax), 2*ndim) # Alternate Lmin, Lmax for each axis
        #aranges = tmp.copy()
        #aranges.resize(6) # Needs to have 6 entries regardless of ndim.
        #print " * Axis ranges (orthorhomic): ",aranges
        #griddata = np.reshape(allcols[ndim], griddim)
        #return griddata, griddim, ndim, aranges

        # DONE

    def readFields(self, filepath, verbose=False):
        """ Wrapper function to parse field grid data from PolyFTS formatted or unformatted files"""
        # First check whether the file is gzipped
        self._testgzipFile(filepath)
        if verbose:
            if self._gzipped:
                print " * File {0} is gzipped".format(filepath)
            else:
                print " * File {0} is not gzipped".format(filepath)
        # Next check whether file is binary formatted PolyFTS file
        self._testBinFile(filepath)
        if verbose:
            if self._BinFile:
                print " * File {0} is an unformatted field file with version #{1}".format(filepath,self.version)
            else:
                print " * File {0} is not an unformatted field file".format(filepath)
        if self._BinFile:
            # Candidate for unformatted PolyFTS data file
            self._readBinFile(filepath)
        else:
            # Else assume that the file is a formatted PolyFTS file
            self._readDatFile(filepath)

        if verbose:
            print " * Format version = {0}".format(self.version)
            print " * Number of fields = {0}".format(self.nfields)
            print " * Number of spatial dimensions = {0}".format(self.Dim)
            if self.Dim == 3:
                print " * Grid dimension = {0}x{1}x{2}".format(*self.griddim)
            elif self.Dim == 2:
                print " * Grid dimension = {0}x{1}".format(*self.griddim)
            elif self.Dim == 1:
                print " * Grid dimension = {0}".format(*self.griddim)
            print " * Total # plane waves = {0}".format(self.NPW)
            print " * K-space data? ",self.kspacedata
            print " * Complex data? ",self.complexdata
