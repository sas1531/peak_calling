import sys,os
from pandas import Series,DataFrame
import pandas as pd
import numpy as np
import re
from datetime import datetime
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import openpyxl
from matplotlib.pyplot import cm
from scipy import stats
pd.options.mode.chained_assignment = None  # default='warn'

########## Part 1 (75 points)

# Please describe the lines indicated below with a "###" symbol.
# Provide your comments in each section immediately following the line or code block
# denoted with "Insert Comments"
# For one line a sentence or two is sufficient, but longer code blocks
# should have 3-5 sentences of explanation for full credit

plotlist=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16] ### Make list of 16 numbers (p/ots)

## Insert Comments ##

## This is a list of all of the output plots
## Each number correpsonds to a type of plot that's created below.


dir='MYC'
if len(sys.argv)>1:
	dir=sys.argv[1]
width=100					### Define width (window size)
if len(sys.argv)>2:			### If you have more than 2 arguments..
	width=int(sys.argv[2]) 	### 	then your width is the 2nd arguemnt

## Insert Comments ##

## Width is the size of the window.
## If the number of arguments that will be passed through the script is
## greater than 2, the second argument is width as an integer.
## This sets up your parameters for the command line before you run the script.


outdir='results-'+dir+'-rmsd-'+str(width)
if not os.path.isdir(outdir):
	os.mkdir(outdir)
dfs={}
peaklist={}
rmsd_th=1.3
rmsd_average=1
height_th=1.0
for filename in os.listdir(dir):
	if filename.endswith('.txt'):
		print(filename)
		count=1
		dfs[filename]=pd.read_table(dir+'/'+filename,header=None) ### Create a dictionary of the files (key) to tables named rmsd+filename (values)

		## Insert Comments ##

		## This chunk creates a dictionary with all of the files and associated data frames.
		## Each filename that ends with .txt is a key in the dictonary dfs.
		## The value assigned to the key is the data in the file that has been read
		## into a data frame using pandas read_table function.
		## The first part of the function is the path and there is no header in the file.


		# Plot 1
		dfs[filename].columns=['pos','intensity']
		if count in plotlist:
			print(count)
			fig,ax = plt.subplots(1,1,figsize=(16,6), sharex=True, sharey=False) 						### Define the amount of subplots (1) and the size, share the x axis (not the y) with all subplots
			ax.set_title(filename[:-4],fontsize=20, fontweight='bold') 									### Label the plot with the cell line name, define font size and weight
			ax.set_xlabel('position',fontsize=20, fontweight='bold') 									### Label the x axis with position, define font size and weight
			ax.set_ylabel('intensity',fontsize=20, fontweight='bold') 									### Label the y axis with intensity, define font size and weight
			ax.plot(dfs[filename]['pos'],dfs[filename]['intensity'],color='darkred',lw=2,alpha=1) 		### Plot points (positoin, intesnity) in dark red, set line width and blending (alpha)
			ax.set_xlim([0,np.max(dfs[filename]['pos'])]) 												### Set the x limit to the maximum position value
			ax.set_ylim([np.min(dfs[filename]['intensity']),1.02*np.max(dfs[filename]['intensity'])]) 	### Set the y limit as 1.02 times the maximum intesntiy value
			fig.savefig(outdir+'/'+filename[:-4]+'.'+str(count)+'.'+'.png',dpi=72,bbox_inches='tight') 	### Save the figure with resolution 72 and fit the entire figure in this box (tight)
			plt.close(fig)

			## Insert Comments ##

			## This creates figure 1 and one subplot based on the file columns posiiton and intensity.
			## The title of the plot is the first component of the filename.
			## The x axis is position and the y asix is intensity (specify labels font size and weight).
			## Plot the position against intensity.
			## Set the x axis limit to the maximum position value
			## Set the y axis limist to 1.02 times that intensity value
			## Save this figure as file name and count (1) as png, with resolution 72 and bbox_inches
			## 'tight' to make sure all plot components are saved (nothing is cut off).


		count=2
		window='flat'
		w=np.ones(2*width+1, 'd') ###
		dfs[filename]['avg']=np.convolve(w/w.sum(),dfs[filename]['intensity'],'same') 	### Calculate the moving average of intensity
		dfs[filename]['sd']=(dfs[filename]['intensity']-dfs[filename]['avg'])**2 		### Calcualte the squared deviation of intensity
		dfs[filename]['msd']=np.convolve(w/w.sum(),dfs[filename]['sd'],'same')			### Calculate the mean squared deviation of intensity
		dfs[filename]['rmsd']=dfs[filename]['msd']**0.5 								### Calculate the root mean squared deviaiton of intensity

		## Insert Comments ##

		## This code finds the root mean square deviation of the data (smooths data):
		## Create an array with specified length filled with ones to be used in the rmsd formula.
		## Divide the sum of intensity from 1 (w) and calculate the moving average of the intensity
		## using convolution (giving you the smoothed values).
		## Calculate the squared deviation from the average.
		## Divide the sum of the squared deviation from 1(w) and calculate the moving average of the mean
		## of the squared deviation using convolution (giving you the smoothed values again).
		## Take the root of the mean squared devitaion to get the root mean square deviation.
		## Record all of the above values under the file name in the dictionary.

		# Plot 2
		if count in plotlist:
			print(count)
			fig,ax = plt.subplots(1,1,figsize=(16,6), sharex=True, sharey=False)
			ax.set_title(filename[:-4],fontsize=20, fontweight='bold')
			ax.set_xlabel('position',fontsize=20, fontweight='bold')
			ax.set_ylabel('rmsd',fontsize=20, fontweight='bold')
			ax.plot(dfs[filename]['pos'],dfs[filename]['rmsd'],color='black',lw=2,alpha=1)
			ax.set_xlim([0,np.max(dfs[filename]['pos'])])
			ax.set_ylim([0,1.02*np.max(dfs[filename]['rmsd'])])
			fig.savefig(outdir+'/'+filename[:-4]+'.'+str(count)+'.'+'-'+str(width)+'-rmsd.png',dpi=72,bbox_inches='tight')
			plt.close(fig)
		count+=1

		# Plot 3
		fig,ax = plt.subplots(1,1,figsize=(6,6), sharex=True, sharey=False)
		ax.set_title(filename[:-4],fontsize=20, fontweight='bold')
		ax.set_xlabel('Intensity',fontsize=20, fontweight='bold')
		n, bins, patches = ax.hist(dfs[filename]['intensity'],color='black',bins=100,lw=2,alpha=1,histtype='step')
		#ax.set_xlim([0,np.max(bins)])
		if count in plotlist:
			print(count)
			fig.savefig(outdir+'/'+filename[:-4]+'.'+str(count)+'.'+'-'+str(width)+'-int-hist.png',dpi=72,bbox_inches='tight')
		count+=1
		max_val=np.max(n)
		th=0.2
		th_=0.1
		th_i=0
		max_i=len(n)
		done=0
		for i in range(len(n)):				### Iterate through the length of list n
			if done==0:						### Set done (can be used to exit loop later)
				if max_val==n[i]:			### If n[i] equals your pre-set max value...
					max_i=i 				### 	set your max value bin (max_i) to i
				if i>max_i: 				### If i is larger than your max_i (bin)...
					if max_val*th_>n[i]:   	### 	AND if n[i] is less than 10% of your max value..
						done=1 				### 	set done to 1 and break out of the loop
				if max_val*th<n[i]: 		### If n[i] is larger than 20% of the max value..
					th_i=i+1 				###		save this bin number (plus one) as your threshold bin


		## Insert Comments ##

		## This finds the threshold of the signal data.

		## Iterate through n:
		## Find the maximum bin count and store which bin it is as max_i (max bin #).
		## If there is a new bin # that is greater than max_i but has less than 10% of the
		## maximum count then break out of the loop (this skips low count bins located after the threshold).
		## If the bin count is greater than 20% of the maximum count then the threhold
		## is that bin # plus 1 (this finds the actual threshold bin = where to draw the line).

		## n is a list of the number of counts in each bin of the histogram
		## max_val = maximum count (bin with the most values)
		## th = threshold
		## th_i = threhold bin #
		## max_i = maximum count bin #


		# Plot 4
		ax.plot([bins[th_i],bins[th_i]],[0,max_val],color='darkred',lw=1,alpha=1)
		# ^this just plots a line at bin[thi_] from the x axis to the max y value...
		if count in plotlist:
			print(count)
			fig.savefig(outdir+'/'+filename[:-4]+'.'+str(count)+'.'+'-'+str(width)+'-int-hist_.png',dpi=72,bbox_inches='tight')
		count+=1
		plt.close(fig)
		intensity_level=bins[th_i]

		# Plot 5
		fig,ax = plt.subplots(1,1,figsize=(6,6), sharex=True, sharey=False)
		ax.set_title(filename[:-4],fontsize=20, fontweight='bold')
		ax.set_xlabel('rmsd',fontsize=20, fontweight='bold')
		n, bins, patches = ax.hist(dfs[filename]['rmsd'],color='black',bins=100,lw=2,alpha=1,histtype='step')
		ax.set_xlim([0,np.max(bins)])
		if count in plotlist:
			print(count)
			fig.savefig(outdir+'/'+filename[:-4]+'.'+str(count)+'.'+'-'+str(width)+'-rmsd-hist.png',dpi=72,bbox_inches='tight')
		count+=1
		max_val=np.max(n)
		th=0.2 ### Threshold value (20% of the maximum)

		## Insert Comments ##

		## This is the intensity threshold percent of the maximum (20%).
		## Only look at regions above this threshold in later analysis.


		th_=0.1
		th_i=0
		max_i=len(n)
		done=0
		for i in range(len(n)): 			### Iterate through the length of list n
			if done==0:						### Set done (can be used to exit loop later)
				if max_val==n[i]:			### If n[i] equals your pre-set max value...
					max_i=i 				### 	set your max value bin (max_i) to i
				if i>max_i: 				### If i is larger than your max_i (bin)...
					if max_val*th_>n[i]: 	### 	AND if n[i] is less than 10% of your max value..
						done=1 				### 	set done to 1 and break out of the loop
				if max_val*th<n[i]: 		### If n[i] is larger than 20% of the max value..
					th_i=i+1 				### 	save this bin number (plus one) as your threshold bin
		rmsd_level=bins[th_i] 				### rmsd_level defined as the threshold (red line location) for the rmsd plot
		#print(rmsd_level)

		## Insert Comments ##

		## Same iteration as previous chunk, used to find the threshold of the signal data.
		## Iterate through n:
		## Find the maximum bin count and store which bin it is as max_i (max bin #).
		## If there is a new bin # that is greater than max_i but has less than 10% of the
		## maximum count then break out of the loop (this skips low count bins located after the threshold).
		## If the bin count is greater than 20% of the maximum count then the threhold
		## is that bin # plus 1 (this finds the actual threshold bin = where to draw the line).


		if rmsd_average==1:
			print('rmsd avg')
			rmsd_mean=np.mean(dfs[filename]['rmsd'][dfs[filename]['rmsd']<=rmsd_level]) 	### Find the mean of rmsd
			rmsd_std=np.std(dfs[filename]['rmsd'][dfs[filename]['rmsd']<=rmsd_level])		### Find the standard deviation of rmsd
			rmsd_level=rmsd_mean+rmsd_std*(-2*np.log(th))**0.5 								### Re-define the rmsd level using the mean and standard deviation

			## Insert Comments ##

			## This loop will always run because rmds_average is previously set to 1:
			## This finds the mean and the standard deviation of rmsd that's less than or
			## equal to the root mean square deviation threshold.
			## The standard deviation and mean are then used to find the new rmsd threshold
			## for the rmsd histogram that is plotted below.

		# Plot 6
		ax.plot([rmsd_level,rmsd_level],[0,max_val],color='darkred',lw=1,alpha=1)
		if count in plotlist:
			print(count)
			fig.savefig(outdir+'/'+filename[:-4]+'.'+str(count)+'.'+'-'+str(width)+'-rmsd-hist_.png',dpi=72,bbox_inches='tight')

		# Plot 7
		count+=1
		ax.plot([rmsd_level*rmsd_th,rmsd_level*rmsd_th],[0,max_val],color='darkred',lw=1,ls='--',alpha=1)
		if count in plotlist:
			print(count)
			fig.savefig(outdir+'/'+filename[:-4]+'.'+str(count)+'.'+'-'+str(width)+'-rmsd-hist__.png',dpi=72,bbox_inches='tight')

		# Plot 8
		count+=1
		plt.close(fig)

		if count in plotlist:
			print(count)
			fig,ax = plt.subplots(1,1,figsize=(16,6), sharex=True, sharey=False)
			ax.set_title(filename[:-4],fontsize=20, fontweight='bold')
			ax.set_xlabel('position',fontsize=20, fontweight='bold')
			ax.set_ylabel('intensity',fontsize=20, fontweight='bold')
			ax.plot(dfs[filename]['pos'],dfs[filename]['intensity'],color='black',lw=2,alpha=1)
			for idx in dfs[filename]['pos'].index:
				if dfs[filename]['rmsd'].ix[idx]>rmsd_level*rmsd_th:
					ax.scatter(dfs[filename]['pos'].ix[idx],1.01*np.max(dfs[filename]['intensity']),color='darkred',s=10,alpha=1)
			ax.set_xlim([0,np.max(dfs[filename]['pos'])])
			ax.set_ylim([np.min(dfs[filename]['intensity']),1.02*np.max(dfs[filename]['intensity'])])
			fig.savefig(outdir+'/'+filename[:-4]+'.'+str(count)+'.'+'-'+str(width)+'_.png',dpi=72,bbox_inches='tight')
			plt.close(fig)
		count+=1

		width2=width
		window2='hanning'
		intensity_max=0
		peak_width_threshold=0.25
		if window2 in ['hanning', 'hamming', 'bartlett', 'blackman']:
			w2=eval('np.'+window2+'(2*width2+1)')
		else:
			window2='flat'
			w2=np.ones(2*width2+1,'d')
		dfs[filename]['int-smooth']=np.convolve(w2/w2.sum(),dfs[filename]['intensity'],'same') ### Calculate teh moving average of intensity and smooth via hanning

		## Insert Comments ##

		## Hanning window is a tapering function that uses weighted cosine (one period of cosine plus 1)
		## - this will be used to smooth the data.
		## Divide the sum of the intensity from the hanning window (w2) and calculate the moving average
		## of intensity using convolution (giving you the smoothed values).
		## Assign these values to the file as a series to be graphed below.


		# Plot 9
		if count in plotlist:
			print(count)
			fig,ax = plt.subplots(1,1,figsize=(16,6), sharex=True, sharey=False)
			ax.set_title(filename[:-4],fontsize=20, fontweight='bold')
			ax.set_xlabel('position',fontsize=20, fontweight='bold')
			ax.set_ylabel('intensity',fontsize=20, fontweight='bold')
			ax.plot(dfs[filename]['pos'],dfs[filename]['int-smooth'],color='black',lw=2,alpha=1)
			ax.set_xlim([0,np.max(dfs[filename]['pos'])])
			ax.set_ylim([np.min(dfs[filename]['intensity']),1.02*np.max(dfs[filename]['intensity'])])
			fig.savefig(outdir+'/'+filename[:-4]+'.'+str(count)+'.'+'-'+window+str(width)+'-int-smooth.png',dpi=72,bbox_inches='tight')
			plt.close(fig)
		count+=1
		dfs[filename]['int-diff']=dfs[filename]['int-smooth'].diff() ###

		## Insert Comments ##

		## The .diff function calculates the difference between two neighboring elements.
		## This is used to find the derivate of the smoothed data (which can then be used to
		## find the maxima and minima).


		# Plot 10
		if count in plotlist:
			print(count)
			fig,ax = plt.subplots(1,1,figsize=(16,6), sharex=True, sharey=False)
			ax.set_title(filename[:-4],fontsize=20, fontweight='bold')
			ax.set_xlabel('position',fontsize=20, fontweight='bold')
			ax.set_ylabel('intensity',fontsize=20, fontweight='bold')
			ax.plot(dfs[filename]['pos'],dfs[filename]['int-diff'],color='black',lw=2,alpha=1)
			ax.set_xlim([0,np.max(dfs[filename]['pos'])])
			fig.savefig(outdir+'/'+filename[:-4]+'.'+str(count)+'.'+'-'+window+str(width)+'-int-diff.png',dpi=72,bbox_inches='tight')
			plt.close(fig)
		count+=1
		dfs[filename]['int-diff-bin']=0 													### Assign all positions in this series as 0
		dfs[filename]['int-diff-bin'][dfs[filename]['int-diff']>0]=1 						### If the derivative is greater than 0, assign the position a 1
		dfs[filename]['int-diff-bin-diff']=dfs[filename]['int-diff-bin'].diff() 			### Calculate the second derivative of each position
		dfs[filename]['int-diff-bin-diff-max']=dfs[filename]['int-diff-bin-diff']			### Create a max second derivative bin
		dfs[filename]['int-diff-bin-diff-min']=dfs[filename]['int-diff-bin-diff']			### Create a min secnod derivative bin
		dfs[filename]['int-diff-bin-diff-max'][dfs[filename]['int-diff-bin-diff']>0]=0 		### If the position is positive, assign it as 0
		dfs[filename]['int-diff-bin-diff-min'][dfs[filename]['int-diff-bin-diff']<0]=0 		### If the position is negative, assign it as 0

		## Insert Comments ##

		## This chunk creates and stores the maximas and minimas of the plot using the derivative values.

		## Create a series and assign 0 to each derivative position.
		## If the derivative value at that position is greater than zero assign it a 1 (only the zero values can
		## be maxima or minima).
		## Calculate the second derivative.
		## If the second derivative is positive, assign it a zero (this means all maximas will be -1) and put it
		## in the max series
		## If the second derivative is negative, assign it a zero (this means all minimas will be 1) and put it
		## in the min series.

		## The result is a series showing the positions of the maximas and minimas as 1 (all others are 0).


		# Plot 11
		if count in plotlist:
			print(count)
			fig,ax = plt.subplots(1,1,figsize=(16,6), sharex=True, sharey=False)
			ax.set_title(filename[:-4],fontsize=20, fontweight='bold')
			ax.set_xlabel('position',fontsize=20, fontweight='bold')
			ax.set_ylabel('intensity',fontsize=20, fontweight='bold')
			ax.plot(dfs[filename]['pos'],dfs[filename]['intensity'],color='darkred',lw=2,alpha=1)
			ax.plot(dfs[filename]['pos'],(-dfs[filename]['int-diff-bin-diff-max'])*1.02*(np.max(dfs[filename]['intensity'])-np.min(dfs[filename]['intensity']))+np.min(dfs[filename]['intensity']),color='black',lw=1,alpha=1)
			ax.set_xlim([0,np.max(dfs[filename]['pos'])])
			ax.set_ylim([np.min(dfs[filename]['intensity']),1.02*np.max(dfs[filename]['intensity'])])
			fig.savefig(outdir+'/'+filename[:-4]+'.'+str(count)+'.'+'-'+window+str(width)+'-int-diff-bin-diff-max.png',dpi=72,bbox_inches='tight')
			plt.close(fig)

		# Plot 12
		count+=1
		if count in plotlist:
			print(count)
			fig,ax = plt.subplots(1,1,figsize=(16,6), sharex=True, sharey=False)
			ax.set_title(filename[:-4],fontsize=20, fontweight='bold')
			ax.set_xlabel('position',fontsize=20, fontweight='bold')
			ax.set_ylabel('intensity',fontsize=20, fontweight='bold')
			ax.plot(dfs[filename]['pos'],dfs[filename]['intensity'],color='darkred',lw=2,alpha=1)
			ax.plot(dfs[filename]['pos'],(dfs[filename]['int-diff-bin-diff-min'])*1.02*(np.max(dfs[filename]['intensity'])-np.min(dfs[filename]['intensity']))+np.min(dfs[filename]['intensity']),color='black',lw=1,alpha=1)
			ax.set_xlim([0,np.max(dfs[filename]['pos'])])
			ax.set_ylim([np.min(dfs[filename]['intensity']),1.02*np.max(dfs[filename]['intensity'])])
			fig.savefig(outdir+'/'+filename[:-4]+'.'+str(count)+'.'+'-'+window+str(width)+'-int-diff-bin-diff-min.png',dpi=72,bbox_inches='tight')
			plt.close(fig)
		count+=1

		peaklist[filename]=[]
		for idx in dfs[filename].index:
			idx_=idx
			if dfs[filename]['int-diff-bin-diff'].ix[idx]==-1: 		### Isolate the maxima (-1)
				peak_pos=dfs[filename]['pos'].ix[idx] 				### Define the maxima peak position at this location
				peak_height=dfs[filename]['intensity'].ix[idx] 		### Define the maxima peak height (intensity) at this location
				signal_to_noise=0 									### Set signal noise to 0
				done=0 												### set done to 0 (use in the next nested loops)

				## Insert Comments ##

				## .index returns the lowest index where the elements appears, use this to
				## iterate through all of the derivative values and find the maxima (-1).
				## When a maxima is found find the position and intensity values and
				## store then as peak_pos and peak_height.
				## This establishes and records the top of the maxima/peak (the actual
				## values/coordinates).


				for i in range(width2): # find the left side of peak mean width height (minus)
					if idx-i in dfs[filename].index:
						if dfs[filename]['intensity'].ix[idx-i]<peak_height/2.0:
							done=1
						if dfs[filename]['int-diff-bin-diff'].ix[idx-i]==1:
							done=1
						if done==0:
							if peak_height<dfs[filename]['intensity'].ix[idx-i]:
								peak_height=dfs[filename]['intensity'].ix[idx-i]
								peak_pos=dfs[filename]['pos'].ix[idx-i]
								idx_=idx-i
				done=0
				for i in range(width2): # find the right side of the peak mean with height (plus)
					if idx+i in dfs[filename].index:
						if dfs[filename]['intensity'].ix[idx+i]<peak_height/2.0:
							done=1
						if dfs[filename]['int-diff-bin-diff'].ix[idx+i]==1:
							done=1
						if done==0:
							if peak_height<dfs[filename]['intensity'].ix[idx+i]:
								peak_height=dfs[filename]['intensity'].ix[idx+i]
								peak_pos=dfs[filename]['pos'].ix[idx+i]
								idx_=idx+i

				i_min_minus=0 # find the left side of peak mean width bin (minus)
				done=0
				for i in range(width2*5):
					if idx_-i in dfs[filename].index:
						if dfs[filename]['int-diff-bin-diff'].ix[idx_-i]==1:
							done=1
						if done==0:
							i_min_minus=i
				i_min_plus=0 # find the right side of peak mean width bin (plus)
				done=0
				for i in range(width2*5):
					if idx_+i in dfs[filename].index:
						if dfs[filename]['int-diff-bin-diff'].ix[idx_+i]==1:
							done=1
						if done==0:
							i_min_plus=i
				local_bgr=(dfs[filename]['intensity'].ix[idx_-i_min_minus]+dfs[filename]['intensity'].ix[idx_+i_min_plus])/2.0
				peak_minus=-1
				i_peak_minus=0
				done=0
				for i in range(i_min_minus):
					if idx_-i in dfs[filename].index:
						if (dfs[filename]['intensity'].ix[idx_-i]-local_bgr)<(peak_height-local_bgr)*peak_width_threshold:
							done=1
						if done==0:
							peak_minus=dfs[filename]['pos'].ix[idx_-i]
							i_peak_minus=i
				peak_plus=-1
				i_peak_plus=0
				done=0
				for i in range(i_min_plus):
					if idx_+i in dfs[filename].index:
						if (dfs[filename]['intensity'].ix[idx_+i]-local_bgr)<(peak_height-local_bgr)*peak_width_threshold:
							done=1
						if done==0:
							peak_plus=dfs[filename]['pos'].ix[idx_+i]
							i_peak_plus=i
				avg=dfs[filename]['intensity'].ix[idx_]
				avg_count=1
				for i in range(i_peak_minus):
					if idx_-1-i in dfs[filename].index:
						avg+=dfs[filename]['intensity'].ix[idx_-1-i]
						avg_count+=1
				for i in range(i_peak_plus):
					if idx_+1+i in dfs[filename].index:
						avg+=dfs[filename]['intensity'].ix[idx_+1+i]
						avg_count+=1
				avg/=avg_count
				rmsd=(dfs[filename]['intensity'].ix[idx_]-avg)**2
				rmsd_count=1
				for i in range(i_peak_minus):
					if idx_-1-i in dfs[filename].index:
						rmsd+=(dfs[filename]['intensity'].ix[idx_-1-i]-avg)**2
						rmsd_count+=1
				for i in range(i_peak_plus):
					if idx_+1+i in dfs[filename].index:
						rmsd+=(dfs[filename]['intensity'].ix[idx_+1+i]-avg)**2
						rmsd_count+=1
				rmsd/=rmsd_count
				rmsd=rmsd**0.5
				signal_to_noise=rmsd/rmsd_level
				if 30<peak_pos and peak_pos<6020:
					peaklist[filename].append((peak_pos,peak_minus,peak_plus,peak_plus-peak_minus,peak_height,peak_height-intensity_level,avg,local_bgr,signal_to_noise))
		f=open(outdir+'/'+filename[:-4]+'-'+window+str(width)+'-peaklist.text','w')
		f.write('pos\twhm-\twhm+\tdelta\theight\theight_norm\tavg\tbgr\tsig_to_noise\n')
		for (pos,minus,plus,delta,height,height_norm,avg,bgr,sig_to_noise) in peaklist[filename]:
			f.write(str(pos)+'\t'+str(minus)+'\t'+str(plus)+'\t'+str(delta)+'\t'+str(height)+'\t'+str(height_norm)+'\t'+str(avg)+'\t'+str(bgr)+'\t'+str(sig_to_noise)+'\n')
		f.close()

		# Plot 13
		fig,ax = plt.subplots(1,1,figsize=(16,6), sharex=True, sharey=False)																				### Define the amount of subplots (1) and the size, share the x axis (not the y) with all subplots
		ax.set_title(filename[:-4],fontsize=20, fontweight='bold')																							### Label the plot with the cell line name, define font size and weight
		ax.set_xlabel('position',fontsize=20, fontweight='bold')																							### Label the x axis with position, define font size and weight
		ax.set_ylabel('intensity',fontsize=20, fontweight='bold')																							### Label the y axis with intensity, define font size and weight
		ax.plot(dfs[filename]['pos'],dfs[filename]['intensity'],color='darkred',lw=2,alpha=1)																### Plot points (position, intesnity) in dark red, set line width and blending (alpha)
		for (pos,minus,plus,delta,height,height_norm,avg,bgr,sig_to_noise) in peaklist[filename]:															### For all files in list of plot caluclations...
			if sig_to_noise>=rmsd_th and height_norm>height_th:																								###		If the signal noise is greater than the threshold AND the normalized height is larger than the threshold height
				ax.scatter([pos],[1.01*height],facecolor='blue',edgecolor='darkblue',s=60,lw=2,alpha=1)														###		Put a blue circle at the top of the peaks (1.02 times the height) on the plot
		ax.set_xlim([0,np.max(dfs[filename]['pos'])])																										### Set the x limit to the maximum position value
		ax.set_ylim([np.min(dfs[filename]['intensity']),1.02*np.max(dfs[filename]['intensity'])])															### Set the y limit from the minimum intensity to 1.02 times the maximum intensity
		if count in plotlist:																																### If the count is in the plotlist (1-16)...
			print(count)																																	###		Print the count...
			fig.savefig(outdir+'/'+filename[:-4]+'.'+str(count)+'.'+'-'+window+str(width)+'-peaks.png',dpi=72,bbox_inches='tight')							###		Save the figure with resolution 72 and fit the entire figure in this box (tight)

		# Plot 14
		count+=1																																			### Add one to the count (for identification in the list and plot saving)
		for (pos,minus,plus,delta,height,height_norm,avg,bgr,sig_to_noise) in peaklist[filename]:															### For all files in list of plot caluclations...
			if sig_to_noise>=rmsd_th and height_norm>height_th:																								###		If the signal noise is greater than the threshold AND the normalized height is larger than the threshold height
				ax.plot([minus,plus],[peak_width_threshold*(height-bgr)+bgr,peak_width_threshold*(height-bgr)+bgr],color='darkblue',lw=2,alpha=1)			###		Plot a horizontal line from minus and plus position at the peak width threshold intensity value (adjusted to background noise?)
		if count in plotlist:																																### If the count is in the plotlist (1-16)...
			print(count)																																	### 	Print the count...
			fig.savefig(outdir+'/'+filename[:-4]+'.'+str(count)+'.'+'-'+window+str(width)+'-peaks-widths.png',dpi=72,bbox_inches='tight')					### 	Save the figure with resolution 72 and fit the entire figure in this box (tight)

		# Plot 15
		count+=1																																			### Add one to the count (for identification in the list and plot saving)
		for (pos,minus,plus,delta,height,height_norm,avg,bgr,sig_to_noise) in peaklist[filename]:															### For all files in list of plot caluclations...
			if sig_to_noise>=rmsd_th and height_norm>height_th:																								### 	If the signal noise is greater than the threshold AND the normalized height is larger than the threshold height
				ax.plot([minus,plus],[bgr,bgr],color='gray',lw=1,alpha=1)																					### 	Plot a horizontal line from minus and plus position at the background noise intensity value
		if count in plotlist:																																### If the count is in the plotlist (1-16)...
			print(count)																																	### 	Print the count...
			fig.savefig(outdir+'/'+filename[:-4]+'.'+str(count)+'.'+'-'+window+str(width)+'-peaks-widths-bgr.png',dpi=72,bbox_inches='tight')				### 	Save the figure with resolution 72 and fit the entire figure in this box (tight)
			plt.close(fig)																																	### Close the figure

		# Plot 16
		count+=1																																			### Add one to the count (for identification in the list and plot saving)
		for (pos,minus,plus,delta,height,height_norm,avg,bgr,sig_to_noise) in peaklist[filename]:															### For all files in list of plot caluclations...
			if sig_to_noise>=rmsd_th and height_norm>height_th:																								### 	If the signal noise is greater than the threshold AND the normalized height is larger than the threshold height
				ax.plot([minus,plus],[avg,avg],color='black',lw=1,alpha=1)																					### 	Plot a horizontal line from minus and plus position at the average intensity value
		if count in plotlist:																																### If the count is in the plotlist (1-16)...
			print(count)																																	### 	Print the count...
			fig.savefig(outdir+'/'+filename[:-4]+'.'+str(count)+'.'+'-'+window+str(width)+'-peaks-widths-bgr-avg.png',dpi=72,bbox_inches='tight')			### 	Save the figure with resolution 72 and fit the entire figure in this box (tight)
			plt.close(fig)																																	### Close the figure
		count+=1																																			### Add one to the count (for identification in the list and plot saving)

		## Insert Comments  ##
		## Generally speaking, how does this block of code differ from the code block ##
		## that plots the first figure

		## The first code block sets up the general frequency plot without any annotation,
		## defining the axes, titles, colors, fonts, weights, etc.
		## This code block prints the exact same frequency data; however it utilizes
		## the peaklist information (created in the nested for loops above) to visualize
		## the peaks, the mean width, the background noise, and average. These values are
		## pulled from the peaklist file that was created separately. Only the peaks that are
		## above the threshold will be shown, visualizing the actual peaks against the
		## background noise.



########## Part 2 (25 points)

# Run the above script at any desired window size for the samples in the "myc/" directory
# Submit with your commented code the final annotated plots for the MYC protein in each cell line (3 plots).
#
# In another text file or below comment on the signal of MYC in each of these cell lines (4-7 sentences).
# Would you expect the signals to differ across the three cell lines?
# How might various window sizes thresholds effect peak calling in this setting?
# (Biological significant or coding difference..?)


# MCF-7 (breast cancer)
# K562 (leukemia)
# HepG2 (liver cancer)
# The MYC family consists of regulator genes and proto-oncogenes that code for
# transcription factors. Up-regulation of myc genes is related to cellular growth
# metabolism for cancerous tumors. Though each of these cell lines is related
# to a different cancer, the general up-regulation of myc shouldn't be significantly
# different between each line as it's a hallmark of cancer growth/metabolism. As such,
# I would predict that the signals should be similar for each cell line. That being said,
# each of the signal frequencies is slightly different, with different peaks and significantly
# different backgrounds. This could imply that there are different levels of up-regulation
# depending on cancer type specifically shown with the number of peaks and the level of
# background noise.
# At window size 100, this appears to be supported where all of the promoter and transcription
# binding site regions have peaks and each cell line has a different level of background noise.
# As you increase the window size, you still see peaks at the beginning of the gene
# suggesting up-regulation. If the window size is too small (<=40) you will get more varied
# peaks across each cell line that aren't as consistent and as the window sizes increase you
# lose peaks.
