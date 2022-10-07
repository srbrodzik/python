import os

indir = '/home/storm/relops/soundings/SMN/gifs'
os.chdir(indir)

for file in os.listdir(indir):
    if file.endswith('.png'):
        (fname,ext) = os.path.splitext(file)
        cmd = 'convert '+fname+ext+' '+fname+'.gif'
        os.system(cmd)
        cmd = 'rm '+file
        os.system(cmd)
        
    
