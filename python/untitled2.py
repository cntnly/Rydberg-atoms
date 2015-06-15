from struct import pack,unpack
import wave
import winsound,sys
import time, datetime
import random
import numpy as np
today=datetime.datetime.now()
today.hour
time.sleep(1)

samplelen=100
noise_output = wave.open('noise2.wav', 'w')
noise_output.setparams((1, 2, samplelen, 0, 'NONE', 'not compressed'))

values = []

for i in range(0, samplelen):
        value = 32767
        packed_value = pack('h', value)
        values.append(packed_value)
     #   values.append(packed_value)

value_str = b''.join(values)
noise_output.writeframes(value_str)

noise_output.close()

winsound.PlaySound('noise2.wav', winsound.SND_FILENAME)


spf = wave.open('noise2.wav','r')

#Extract Raw Audio from Wav File
signal = spf.readframes(-1)
signal = np.fromstring(signal, 'Int16')
fs = spf.getframerate()

#If Stereo
if spf.getnchannels() == 2:
    print ('Just mono files')
    sys.exit(0)


Time=np.linspace(0, len(signal)/fs, num=len(signal))

figure(1);clf()
title('Signal Wave...')
plot(Time,signal,'wo')
show()