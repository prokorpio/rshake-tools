# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 12:13:21 2022

@author: User
"""

import time
from gtts import gTTS

from pydub import AudioSegment
from pydub.playback import play

import threading
from tkinter import *
from tkinter import messagebox as mb

def intensityConvert(intensity): #convert string to number for text-to-voice
    if intensity == "I":
        intensityNum = 1
    elif intensity == "II-III":
        intensityNum = 2
    elif intensity == "IV":
        intensityNum = 4
    elif intensity == "V":
        intensityNum = 5
    elif intensity == "VI":
        intensityNum = 6
    elif intensity == "VII":
        intensityNum = 7
    elif intensity == "VIII":
        intensityNum = 8
    elif intensity == "IX":
        intensityNum = 9
    elif intensity == "X+":
        intensityNum = 10
    else:
        intensityNum = 0
        
    return intensityNum

def threadSound(): #thread for playing alert sound
    for _ in range(1):
        play(AudioSegment.from_mp3("alarm.mp3"))
        play(AudioSegment.from_mp3("eq_intensity.mp3"))
        play(AudioSegment.from_mp3("eq_displacement.mp3"))
        time.sleep(1)

def threadDialog(intensity,displacement): # thread for displaying alert box
    root = Tk()
    root.withdraw()
    mb.showwarning('Earthquake Alert',
                   'Alert. Earthquake ongoing at Intensity %s, Displacement is %.2f centimeters' % (intensity,displacement),
                   parent=root)
    root.destroy()

def alarm(intensity, displacement, axis):
    # Convert intensity string into number
    intensityNum = intensityConvert(intensity)
    # create audio file for intensity
    if intensity=="II-III":
        myobj = gTTS(text="Alert. Earthquake ongoing at intensity two or three", lang="en", slow=False)
    else:
        myobj = gTTS(text="Alert. Earthquake ongoing at intensity %.0f" % intensityNum, lang="en", slow=False)
    myobj.save("./eq_intensity.mp3")
    # create audio file for displacement
    myobj = gTTS(text="Displacement is %.2f centimeters along"+axis % displacement, lang="en", slow=False)
    myobj.save("./eq_displacement.mp3")

    # Create separate threads for alert sounds and dialog box
    alertThread1 = threading.Thread(target=threadSound)
    alertThread1.start()
    alertThread2 = threading.Thread(target=threadDialog(intensity,displacement))
    alertThread2.start()

    # alertThread1.join()
    # alertThread2.join()

if __name__ == "__main__":
    intensity = "VII"
    displacement = 4.32
    alarm(intensity, displacement)
