#! /usr/bin/python
import math;
import argparse;
import matplotlib.pyplot as plt
import string;
import re;
import numpy as np;
plt.switch_backend('agg')

parser = argparse.ArgumentParser();
parser.add_argument('filename', type=str);
parser.add_argument('--fixed_angle', type=float, default=-10);
parser.add_argument('--mode', type=str, default='ppi');
parser.add_argument('--saveimg', type=str, default='');
args = parser.parse_args();

elev_rate = 1.5
elev_start = 0

#filename = 'sea20171029_235501_2.txt'
filename = args.filename;

if (args.fixed_angle == -10):
    fnp = re.compile('.+_\d+_(\d)');
    sn = fnp.search(args.filename);
    sweep_number = int(sn.groups()[0]) - 1;
    angle_list = [0.5, 0.8, 1.5, 2.0, 3.0, 3.75, 4.75]
    args.fixed_angle = angle_list[sweep_number];


with open(filename, 'r') as f:
    read_data = f.read()

read_data = read_data.replace('\n','');
#read_data = read_data.replace('#','\n#');

all_lines = string.split(read_data, '#');

pattern = r'\s*(\d+).+\s(\d+:\d+:\d+\.\d+)\s+Az:\s*(\d+\.\d+)\s+El:\s+(-*\d+\.\d+)\s+Pitch:\s+(-*\d+\.\d+)\s+Roll:\s+(-*\d+\.\d+)\s+Head:\s+(\d+\.\d+)\s+Vel:\s+(-*\d+\.\d+)\s+deg/s\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+Tr:\s+(-*\d+\.\d+)\s+El_or:\s+(-*\d+\.\d+)'
prog = re.compile(pattern);

time_pattern = r'(\d+):(\d+):(\d+\.\d+)';
time_prog = re.compile(time_pattern);

ray_index = [];
ray_time_str = [];
ray_time = [];
azimuth = [];
elevation = [];
pitch = [];
roll = [];
heading = [];
az_vel = [];
el_vel = [];
pitch_vel = [];
roll_vel = [];
hdg_vel = [];
ped_az = [];
ped_el = [];
calc_az = [];
calc_el = [];
fixed_angle = [];

def generate_rotation_matrix(pitch, roll, yaw):
    cP = math.cos(math.pi/180.0 * pitch);
    sP = math.sin(math.pi/180.0 * pitch);
    cR = math.cos(math.pi/180.0 * roll);
    sR = math.sin(math.pi/180.0 * roll);
    cT = math.cos(math.pi/180.0 * yaw);
    sT = math.sin(math.pi/180.0 * yaw);
    Ax = np.matrix([\
            [1,  0,  0],\
            [0, cR, -sR],\
            [0, sR,  cR]]);
    Ay = np.matrix([\
            [ cP, 0, sP],\
            [  0, 1,  0],\
            [-sP, 0, cP]]);
    Az = np.matrix([\
            [cT, -sT, 0],\
            [sT,  cT, 0],\
            [ 0,   0, 1]]);
    return Az*Ay*Ax

def generate_unit_vector(az, el):
    cAz = math.cos(math.pi/180.0*az);
    sAz = math.sin(math.pi/180.0*az);
    cEl = math.cos(math.pi/180.0*el);
    sEl = math.sin(math.pi/180.0*el);

    return np.matrix([\
            [cAz*cEl],
            [sAz*cEl],
            [sEl]]);

def unit_vector_to_azel(v):
    az = 180.0/math.pi*(math.atan2(v[1], v[0]));
    el = 180.0/math.pi*math.atan2(v[2],\
            math.sqrt(pow(v[0],2) + pow(v[1],2)));
    if (az < 0):
        az = az + 360;
        
    return (az, el)

rhi_elevation = elev_start;
rhi_last_update_time = 0;

for line in all_lines:
    result = prog.search(line);
    if result != None:
        ray_index.append(int(result.groups()[0]))
        ray_time_str.append(result.groups()[1])
        rts = result.groups()[1];
        time_result = time_prog.search(rts);
        rt = int(time_result.groups()[0]) * 60 * 60 +\
                int(time_result.groups()[1]) * 60 +\
                float(time_result.groups()[2]);
        ray_time.append(rt);
        azimuth.append(float(result.groups()[2]))
        elevation.append(float(result.groups()[3]))
        pitch.append(float(result.groups()[4]))
        roll.append(float(result.groups()[5]))
        heading.append(float(result.groups()[6]))
        az_vel.append(float(result.groups()[7]))
        el_vel.append(float(result.groups()[8]))
        pitch_vel.append(float(result.groups()[9]))
        roll_vel.append(float(result.groups()[10]))
        hdg_vel.append(float(result.groups()[11]))
        ped_az.append(float(result.groups()[12]))
        ped_el.append(float(result.groups()[13]))
        if (rhi_last_update_time == 0):
            rhi_elevation = elev_start;
        else:
            rhi_elevation += (elev_rate) * (rt - rhi_last_update_time);
        rhi_last_update_time = rt;
        if (args.mode == 'rhi'):
            fixed_angle.append(rhi_elevation)
        else:
            fixed_angle.append(args.fixed_angle)
        
        A = generate_rotation_matrix(pitch[-1], roll[-1], heading[-1]);
        PedCartesian = generate_unit_vector(azimuth[-1], fixed_angle[-1]); 
        
        EarthCartesian = A.T * PedCartesian;
        
        c_az, c_el = unit_vector_to_azel(EarthCartesian);
      
        calc_az.append(c_az);
        calc_el.append(c_el);
       
ray_time, azimuth, elevation, pitch, roll, pitch_vel, roll_vel, hdg_vel, ped_az, ped_el, el_vel, calc_az, calc_el, fixed_angle =\
  zip(*sorted(zip(ray_time, azimuth, elevation, pitch, roll, pitch_vel, roll_vel, hdg_vel, ped_az, ped_el, el_vel, calc_az, calc_el, fixed_angle)));

ray_time = ray_time - np.min(ray_time);

if (args.mode == 'rhi'):
	el_label = '$\phi_\mathrm{e}$ (%.1f - %.1f)' % (np.min(elevation), np.max(elevation))
	fixed_angle_label = '$\phi$'
else:
	el_label = r'$\phi_\mathrm{e}$ (%.2f)' % (np.mean(elevation))
	fixed_angle_label = '$\phi$ (%.1f)' % (args.fixed_angle)
	
roll_label = 'Roll (%.1f)' % (np.max(roll) - np.min(roll))
pitch_label = 'Pitch (%.1f)' % (np.max(pitch) - np.min(pitch))

figscale = 1.25;
fig = plt.figure(tight_layout=True, figsize=(4.5*figscale,5.75*figscale), dpi=100)
plt.subplot(311)

plt.plot(ray_time, pitch, label=pitch_label);
plt.plot(ray_time, roll, label=roll_label);
#plt.plot(ray_time, np.abs(pitch_vel), label='pitch_vel');
#plt.plot(ray_time, np.abs(roll_vel), label='roll_vel');

plt.title(filename);
plt.xlabel('Time (s)');
plt.ylabel('Angle (deg)');
plt.grid();
plt.legend(fontsize='small', fancybox=True, loc='best', ncol=2);
angles = np.concatenate((roll, pitch));
plt.axis([0, max(ray_time), min(angles) * 1.25, max(angles) * 1.25])

plt.subplot(312)
plt.plot(ray_time, elevation, label=el_label);
plt.plot(ray_time, ped_el, label='$\phi_\mathrm{p}$');
#plt.plot(ray_time, el_vel, label='el_vel');
plt.plot(ray_time, fixed_angle, label=fixed_angle_label, linestyle='dashed');
plt.plot(ray_time, calc_el, label='$\phi_\mathrm{cmd}$');

plt.xlabel('Time (s)');
plt.ylabel('Angle (deg)');
plt.grid();
plt.legend(fontsize='x-small', fancybox=True, loc='best', ncol=2);
angles = np.concatenate((ped_el, elevation));
midpoint = 0.5 * (max(angles) + min(angles))
span = (max(angles) - min(angles)) * 0.5 * 1.25
plt.axis([0, max(ray_time), midpoint - span, midpoint + span])

plt.subplot(313)
#plt.plot(ray_time, azimuth, label='azimuth');
#plt.plot(ray_time, calc_az, label='calc_az');
#plt.plot(ray_time, ped_az, label='ped_az');
error1 = np.array(elevation) - np.array(fixed_angle);
error2 = np.array(ped_el) - np.array(calc_el);
e1_label = '$\phi_\mathrm{e}-\phi (\sigma=%.2f)$'%(np.std(error1));
e2_label = '$\phi_\mathrm{p}-\phi_\mathrm{cmd} (\sigma=%.2f)$'%(np.std(error2));
plt.plot(ray_time, error1, label=e1_label);
plt.plot(ray_time, error2, label=e2_label);
plt.xlabel('Time (s)');
plt.ylabel('Error (deg)');
plt.grid();
plt.legend(fontsize='small', fancybox=True, loc='best', ncol=2);
plt.axis([0, max(ray_time), -0.5, 0.5]);

if (len(args.saveimg) > 0):
	plt.savefig(args.saveimg);
else:	
	plt.show();


