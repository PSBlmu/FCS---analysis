# -*- coding: utf-8 -*-
"""
Created on Sat Feb 02 09:35:25 2019

@author: LAdmin
"""

# Background (manually)
import numpy as np
import struct

### Functions to read .sm-files ###

def _get_sm_offset(file_name):
    """
    Not all sm files are exactly the same apperently. Why would you define a standard? Standards are boring.
    It's better to let everybody waste a couple of days on reverse engineering the format.
    p.s. this function returns the number of bytes by which the data is shift from our own microscopes, aka the standard.
    """
    with open(file_name, 'rb') as file_handle:
        file_handle.seek(0xa5)
        char = struct.unpack('>c', file_handle.read(struct.calcsize('>c')))[0]
        if char == '2':
            return 1
        else:
            return 0

def _get_number_photons(file_name,fmt='>Qxxx?'):
    # After the text 'Arrival Time Counter' is an unsigned long that gives us the number of characters in the run
    with open(file_name, 'rb') as file_handle:
        file_handle.seek(0x2e)
        fmt_no_char = '>L'
        size_fmt_no_char = struct.calcsize(fmt_no_char)
        # In the no_char, we have to compensate for the fact that the ulong with the number of char includes part of the header so we substract that.
        no_char = 0x2e + 4 - (0xa5 + _get_sm_offset(file_name)) + struct.unpack(fmt_no_char, file_handle.read(size_fmt_no_char))[0]
    no_photons = no_char // struct.calcsize(fmt)
    return no_photons
        

def _identify_channels(data,period_shift,alex_period,start_green,start_red):
    """ Identifies the real channel for every photon, using knowledge of
    which laser was on, and where the photon was detected.
    Modifies data in place, so no return needed.
    """
    # The arrival time of the photon in the bin
    arrival_time_binned = np.mod(data['ArrivalTime'] - period_shift, alex_period)
    # if self.setting.red_laser_on <= arrival_time_binned < self.setting.red_laser_off:
    #     active_laser = False
    # else:
    #     active_laser = True
    print(np.min(arrival_time_binned), np.max(arrival_time_binned))
    active_laser = np.logical_and(
        arrival_time_binned >= start_green,
        arrival_time_binned < start_red,
    )
    # active_laser = arrival_time_binned < self.setting.first_half  # if true, green laser active, else red
    # Replace FinalChannel by a value from 0-3, where 0 -> DD, 1 -> AA, 2 -> DA, 3 -> AD
    data['FinalChannel'][
        np.logical_and(data['Channel'], active_laser)] = 0  # Laser green, detection also happend in green
    data['FinalChannel'][
        np.logical_and(~data['Channel'], ~active_laser)] = 1  # Laser red, detection also happend in red
    data['FinalChannel'][
        np.logical_and(~data['Channel'], active_laser)] = 2  # Laser green, detection happend in red
    data['FinalChannel'][
        np.logical_and(data['Channel'], ~active_laser)] = 3  # Laser red, detection happend in green


def read_to_np(file_name, time_resolution, channel_flip,period_shift,alex_period,start_green,start_red):
    """
    Reads an sm file and returns a structured array with all the info.
    The array contains:
    'ArrivalTime': arrival time of the photon in seconds
    'FinalChannel': 0 -> DD, 1 -> AA, 2 -> DA, 3 -> AD
    'Channel': bool, if true, green laser active, else red
    """
    # todo: Does not yet support kinky .sm files
    record_dtype = np.dtype([
        ('ArrivalTime', '>u8'),
        ('FinalChannel', '>u2'),
        # Risky proposition. We assign the paddingbytes to this channel, so we get the right columns in our array. It does not exist in the file itself! See smfile documentation.
        ('Channel', '>u2'),
    ])
    with open(file_name, 'rb') as file_handle:
        no_photons = _get_number_photons(file_name)
        file_handle.seek(0xa5 + _get_sm_offset(file_name))
        # Instead of python struct, read to array using numpy fromfile.
        data = np.fromfile(file_handle, dtype=record_dtype, count=no_photons)

    new_data = np.empty(data.shape,
                        dtype=[('ArrivalTime', 'float64'), ('FinalChannel', 'int16'), ('Channel', 'bool')])
    new_data['ArrivalTime'] = data['ArrivalTime'].astype('float64', copy=False) * time_resolution
    new_data['FinalChannel'] = data['FinalChannel']
    new_data['Channel'] = data['Channel'] != 1 if channel_flip else data['Channel']
    _identify_channels(new_data,period_shift,alex_period,start_green,start_red)
    return new_data
    
### Functions to calculate background ###
    

def get_binned_data(photons,binning_time):
    channel_photons = _channels_separate(photons)
    binsize = binning_time
    time_data = {}
    time_dict = {}
    print(channel_photons)
    for key in channel_photons:
        max = np.ceil(channel_photons[key].max() / binsize) * binsize
        time = np.arange(0, max, binsize)
        #NOTE: np.digitize is MUCH slower then np.searchsorted!!
        #todo seems to work this binning but needs checking / confirmation
        inds = np.searchsorted(time, channel_photons[key]) - 1

        #todo fix dis hack where photon arrival time 0.0 ends up in bin -1
        time_data[key] = np.bincount(inds[inds >= 0])

        #combine time array to one?
        time_dict[key] = time
        
    return time_data, time_dict

def get_binned_data_color(photons,binning_time,color):
    channel_photons = _channels_separate_color(photons,color)


    binsize = binning_time
#    time_data = {}
#    time_dict = {}
    
#    for key in channel_photons:
    max = np.ceil(channel_photons.max() / binsize) * binsize
    time = np.arange(0, max, binsize)
    #NOTE: np.digitize is MUCH slower then np.searchsorted!!
    #todo seems to work this binning but needs checking / confirmation
    inds = np.searchsorted(time, channel_photons) - 1

    #todo fix dis hack where photon arrival time 0.0 ends up in bin -1
    time_data = np.bincount(inds[inds >= 0])

    #combine time array to one?
    time_dict = time
        
        
    return time_data, time_dict

#--------------------------------------------------------------------------
# Private API
#--------------------------------------------------------------------------

def _channels_separate(photon_data):
    channel_photons = {}
    _channel = {0: 'DD', 1: 'AA', 2: 'DA', 3: 'AD'}
    
    for ch in np.unique(photon_data['FinalChannel']):
        val = photon_data['ArrivalTime'][photon_data['FinalChannel'] == ch]
        key = _channel[ch]
        channel_photons[key] = val

    return channel_photons

def _channels_separate_color(photon_data, color):
    channel_photons = {}
    _channel = {0: 'DD', 1: 'AA', 2: 'DA', 3: 'AD'}
    
    for ch in np.unique(photon_data['FinalChannel']):
        val = photon_data['ArrivalTime'][photon_data['FinalChannel'] == ch]
        key = _channel[ch]
        channel_photons[key] = val
      
    if color == "green":
        channel_photons = np.concatenate((channel_photons['DD'], channel_photons['AD']), axis=0)
    elif color == "red":
        channel_photons = np.concatenate((channel_photons['DA'], channel_photons['AA']), axis=0)
    else:
        channel_photons = -1

    return channel_photons
    
def script_bg(photon_data,threshold,binning_time):
    # Fuck traits.
#    if isinstance(threshold, _Undefined):
#        return

#    photons = self.data_manager.datasets[0].get_photons()
    photons = photon_data

#    bg_data = BGData(photons, params=self.settings)
    time_data, time = get_binned_data(photons,binning_time)
    
    keys = ['DD', 'DA', 'AA']

    backgrounds = []
    for key in keys:
        data = time_data[key]
        bg_mean = np.mean(data[data < threshold])
#        setattr(self.corr, 'bkg_{}'.format(key), bg_mean)
        backgrounds.append(bg_mean)
    
    return backgrounds