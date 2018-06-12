import os
import numpy
from matplotlib import pyplot as plt, lines
from scipy import interpolate
from math import *
import datetime
#import dictionaryV1


_SAMPLENAME = 'S1'
_SAMPLEDICT = {
	'condTest' : {
		'zoneTest' : {
		#	#'Bkgd' : [145823],
		#	'VB' : [146139],
		#	'Ti' : [146140],
		#	'TiBkgd' : [146141],
			'Cu' : [146144, 146145, 146148],
		#	'CuBkgd' : [146149],
		#	'Au' : [146026]#,
		#	'AuBkgd' : [146153]
		#},
		#'ZC5' : {
			#'VB' : [146154],
			#'Ti' : [146155],
			#'TiBkgd' : [146156],
			#'Cu' : [146157, 146161, 146162],
			#'CuBkgd' : [146149]
		#},
		#'Clean' : {
		#	'VB' : [146164],
		#	'Ti' : [146165],
		#	'TiBkgd' : [146166]
		}
	}
}


#test


def loadDATA(NPZlist):
	# Sums NPZ arrays from the NPZlist and normalizes by the count time throughout the list
	global _DATA
	
	_DATA = numpy.zeros((len(NPZlist), 1024, 1024))
	for NPZnumber in NPZlist:
		ImPath = os.path.realpath("C:\\Users\\Axl\\Desktop\\testi06\\" + str(NPZnumber) +'_PCOImage' + '\\'+ str(NPZnumber) + '_PCOImage.npz')
		print(str(NPZnumber) + '_PCOImage.npz loaded')
		_DATA += numpy.load(os.path.join(ImPath))["data"]
	_DATA = _DATA.sum(0)/_PEAKCOUNTTIME

def plotDATA():
	# Used for debugging
	global _DATA, _PROFILEUP, _PROFILEDOWN, _PEAKPROFILE
	
	line1 = _PROFILEUP
	line2 = _PROFILEDOWN
	fig, (ax1, ax2) = plt.subplots(1,2)
	lineUp = lines.Line2D(line1[0], line1[1], lw=2, color='red', axes=ax1)
	lineDown = lines.Line2D(line2[0], line2[1], lw=2, color='red', axes=ax1) #[X1, X2], [Y1, Y2]
	ax1.add_line(lineUp)
	ax1.add_line(lineDown)
	im = ax1.imshow(_DATA)
	ax2.plot(_PEAKPROFILE[:,1], 'r')
	plt.show()#block=False)
	#plt.pause(1)
	#plt.close()

def GetProfile():
	# exptrait un numpy.ndarray du profile ainsi que les deux profiles extremes pour les tracer si besoin 
	global _DATA, _PEAKPROFILE, _PCONFIG, _PROFILEUP, _PROFILEDOWN
	
	xx = numpy.arange(1024)
	yy = numpy.arange(1024)
	_F = interpolate.RectBivariateSpline(xx, yy, _DATA)
	ndata = int(sqrt((_PCONFIG[0]-_PCONFIG[1])**2+(_PCONFIG[2]-_PCONFIG[3])**2))
	_PEAKPROFILE = numpy.column_stack((numpy.linspace(0,ndata-1,ndata),numpy.zeros(ndata)))
	print(_PEAKPROFILE.shape)
	dx = (_PCONFIG[3]-_PCONFIG[2])/sqrt((_PCONFIG[3]-_PCONFIG[2])**2+(_PCONFIG[1]-_PCONFIG[0])**2)
	dy = (_PCONFIG[1]-_PCONFIG[0])/sqrt((_PCONFIG[3]-_PCONFIG[2])**2+(_PCONFIG[1]-_PCONFIG[0])**2)
	ii = 0
	for di in numpy.arange(-_PCONFIG[4], _PCONFIG[4]+1):
		xx = numpy.linspace(_PCONFIG[0]+di*dx, _PCONFIG[1]+di*dx, ndata)
		yy = numpy.linspace(_PCONFIG[2]-di*dy, _PCONFIG[3]-di*dy, ndata)
		if ii == 0:
			_PROFILEUP = [_PCONFIG[0]+di*dx,_PCONFIG[1]+di*dx], [_PCONFIG[2]-di*dy, _PCONFIG[3]-di*dy]
		_PROFILEDOWN = [_PCONFIG[0]+di*dx,_PCONFIG[1]+di*dx], [_PCONFIG[2]-di*dy, _PCONFIG[3]-di*dy]
		_PEAKPROFILE[:,1] += _F.ev(yy,xx)
		ii += 1

def plotPeakProfile():
	# Used for debugging
	
	
	fig, ax1 = plt.subplots()
	ax1.plot(_PEAKPROFILE[:,1], 'r')
	plt.show()#block=False)
	#plt.pause(1)
	#plt.close()

def fetchDotDat(NPZlist):

	global _PEAKCOUNTTIME, _PEAKDOTDAT, _ROOTPATH, _FLAGSAVE, _KPEAK, anaPath, Kcondition
	
	_PEAKCOUNTTIME = 0.
	for NPZnumber in NPZlist:
		DatPath = os.path.realpath(_ROOTPATH + '\\' + str(NPZnumber) +'.dat' + '\\')
		if os.path.isfile(DatPath):	
			datFile = open(DatPath, "r")
			lines=datFile.readlines()
			STV = int(lines[35].split("=")[1].split(".")[0])
			#posX = lines[49].split("=")[1]
			#posY = lines[50].split("=")[1]
			if (_KPEAK == 'VB' and STV == 198) or (_KPEAK == 'Au' and STV == 113) or (_KPEAK == 'Cu' and STV == 123) or (_KPEAK == 'Ti' and STV == 163) or (_KPEAK == 'AuBkgd' and STV == 103) or (_KPEAK == 'CuBkgd' and STV == 133) or (_KPEAK == 'TiBkgd' and STV == 153) == True:
				for line in lines[56:]:
					posX = line.split('\t')[0]
					posY = line.split('\t')[1]
					Ppeem = line.split('\t')[4]
					PO2 = line.split('\t')[5]
					countTime = line.split('\t')[8]
					_PEAKCOUNTTIME += float(countTime)
					if _FLAGSAVE == True:
						if Kcondition == 'init':
							pass
						else:
							_PEAKDOTDAT.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(NPZnumber, STV, posX, posY, Ppeem, PO2, countTime))
			else:
				if _FLAGSAVE == True:
					error = open(anaPath + "\\errors.txt", 'a+')
					error.write('{} : STV is {} for {}\n'.format(NPZnumber, STV, Kpeak))
								

if __name__ == "__main__":
	
	#_SAMPLENAME, _SAMPLEDICT, = dictionaryV1.sampleC()
	
	_XSTART = 286
	_XEND = 773
	_YSTART = 662
	_YEND = 156
	_WIDTH = 10
	_PCONFIG = [_XSTART, _XEND, _YSTART, _YEND, _WIDTH]
	
	
	_FLAGSAVE = False
	_FLAGPLOT = False
	
	_ROOTPATH = os.path.realpath("C:\\Users\\Axl\\Desktop\\testi06\\")
	anaName = '{}_{}'.format(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M'), _SAMPLENAME)
	anaPath = os.path.realpath(_ROOTPATH + "\\{}".format(anaName))
	
	X1, X2 = 50, 670
	print('X1:', X1-20, X1+20, 'X1:', X2-20, X2+20)
	ndata = int(sqrt((_PCONFIG[0]-_PCONFIG[1])**2+(_PCONFIG[2]-_PCONFIG[3])**2))
	_PEAKPROFILE = numpy.column_stack((numpy.linspace(0,ndata-1,ndata),numpy.zeros(ndata)))
	
	
	Kcondition = 'init'
	
	_KPEAK = 'Au'
	AuBkgd = [145963]
	_AuBkgdPROFILE = numpy.column_stack((numpy.linspace(0,ndata-1,ndata),numpy.zeros(ndata)))
	CuBkgd = [145962]
	_CuBkgdPROFILE = numpy.column_stack((numpy.linspace(0,ndata-1,ndata),numpy.zeros(ndata)))
	
	fetchDotDat(AuBkgd)
	loadDATA(AuBkgd)
	GetProfile()
	Y1, Y2 = _PEAKPROFILE[:,1][X1-15:X1+15].mean(), _PEAKPROFILE[:,1][X2-15:X2+15].mean()
	i = 0
	
	for i in range(len(_PEAKPROFILE[:,1])):
		_PEAKPROFILE[:,1][i] = _PEAKPROFILE[:,1][i]-((Y2-Y1)/(X2-X1))*i
	_PEAKPROFILE[:,1] = _PEAKPROFILE[:,1]*(1000/_PEAKPROFILE[:,1][X1-15:X1+15].mean())
	
	_AuBkgdPROFILE = _PEAKPROFILE
	
	print(_AuBkgdPROFILE)
	
	#plt.legend()
	#plt.show()
	
	_KPEAK = 'Cu'
	fetchDotDat(CuBkgd)
	loadDATA(CuBkgd)
	GetProfile()
	Y1, Y2 = _PEAKPROFILE[:,1][X1-15:X1+15].mean(), _PEAKPROFILE[:,1][X2-15:X2+15].mean()
	i = 0
	for i in range(len(_PEAKPROFILE[:,1])):
		_PEAKPROFILE[:,1][i] = _PEAKPROFILE[:,1][i]-((Y2-Y1)/(X2-X1))*i
	_PEAKPROFILE[:,1] = _PEAKPROFILE[:,1]*(1000/_PEAKPROFILE[:,1][X1-15:X1+15].mean())
	
	_CuBkgdPROFILE = _PEAKPROFILE
	plt.plot(_PEAKPROFILE[:,1], label='extract bkgd')
	plt.plot(_CuBkgdPROFILE[:,1], label='corr bkgd')
	for Kcondition, Dzone in _SAMPLEDICT.items():
		for Kzone, Dpeak in Dzone.items():
			for Kpeak, VnpzList in Dpeak.items():
				print(Kcondition + '\t' + Kzone + '\t' + Kpeak)
				
				_KPEAK = Kpeak			
				peakPath = os.path.realpath(anaPath + '\\' + Kcondition + '\\' + Kzone + '\\' + Kpeak)
				
				if _FLAGSAVE == True:
					try:
						os.makedirs(peakPath)
					except OSError:
						if not os.path.isdir(peakPath):
							raise
					_PEAKDOTDAT = open(os.path.join(peakPath + '\\{}-{}-{}-{}.list'.format(_SAMPLENAME, Kcondition, Kzone, Kpeak)), 'a+')
					_PEAKDOTDAT.write('PCOnumber\tSTV\tposX\tposY\tPpeem\tPO2\tcounttime\n')
				
				fetchDotDat(VnpzList)
				loadDATA(VnpzList)
				GetProfile()
				
				if _FLAGPLOT == True:
					plotDATA()
					#plotPeakProfile()
				
				# normalisation: tilt and offset
				#X1, X2 = 80, 670
				Y1, Y2 = _PEAKPROFILE[:,1][X1-15:X1+15].mean(), _PEAKPROFILE[:,1][X2-15:X2+15].mean()
				i = 0
				for i in range(len(_PEAKPROFILE[:,1])):
					_PEAKPROFILE[:,1][i] = _PEAKPROFILE[:,1][i]-((Y2-Y1)/(X2-X1))*i
				_PEAKPROFILE[:,1] = _PEAKPROFILE[:,1]*(1000/_PEAKPROFILE[:,1][X1-15:X1+15].mean())
				
				
				plt.plot(_PEAKPROFILE[:,1], label = 'data corr')
				plt.plot(_PEAKPROFILE[:,1] - _AuBkgdPROFILE[:,1], label = 'data substract')
				
				#PEAKSUBSTRACT =
				
				if _KPEAK == "Au" :
					_PEAKPROFILE = numpy.column_stack((_PEAKPROFILE, _AuBkgdPROFILE[:,1],_PEAKPROFILE[:,1] - _AuBkgdPROFILE[:,1]))
					#print('yep')
				elif _KPEAK == "Cu" :
					_PEAKPROFILE = numpy.column_stack((_PEAKPROFILE, _CuBkgdPROFILE[:,1],_PEAKPROFILE[:,1] - _CuBkgdPROFILE[:,1]))
					
				plt.legend()
				plt.show()
					
				if _FLAGPLOT == True:
					plotDATA()
					#plotPeakProfile()
					
				if _FLAGSAVE == True:
					_PROFILEDOTDAT = open(os.path.join(peakPath + '\\{}-{}-{}-{}.profile'.format(_SAMPLENAME, Kcondition, Kzone, Kpeak)), 'a+')
					numpy.savetxt(_PROFILEDOTDAT, _PEAKPROFILE, header = 'Pixel\t{}-{}-{}-{}'.format(_SAMPLENAME, Kcondition, Kzone, Kpeak), delimiter = '\t', comments='')#, fmt=['%d', '%f'])

