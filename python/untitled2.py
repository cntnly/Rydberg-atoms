import struck
import wave
import winsound,sys
import time, datetime
today=datetime.datetime.now()
today.hour
time.sleep(1)

samplelen=44100
noise_output = wave.open('noise2.wav', 'w')
noise_output.setparams((2, 2, samplelen, 0, 'NONE', 'not compressed'))

values = []

for i in range(0, samplelen):
        value = int(32767)
        packed_value = pack('h', value)
        values.append(packed_value)
        values.append(packed_value)

value_str = b''.join(values)
noise_output.writeframes(value_str)

noise_output.close()

winsound.PlaySound('noise2.wav', winsound.SND_FILENAME)
