import pylab as pl
from numpy import sin, cos, pi


class RealtimeModel(object):
    def __init__(self, values, Konstant, T0, T1, fCarrier, fAudio):
        self.values = values
        self.Konstant = Konstant
        self.T0 = T0
        self.T1 = T1
        self.fCarrier = fCarrier
        self.fAudio = fAudio

    def SinwaveformGenerator(self, arg):
        Tnext = ((self.Konstant*self.T1)*2)-self.T0
        self.values.append(Tnext)
        self.T0 = self.T1
        self.T1 = Tnext

    def FMwaveformGenerator(self, arg):
        audioInt = -cos(2*pi*self.fAudio*self.T0)
        freqMod = sin(2*pi*self.fCarrier*self.T0 + 2*pi*1*audioInt)
        self.values.append(freqMod)
        self.T0 = self.T1
        self.T1 = freqMod


class RealtimeView(object):
    def __init__(self):
        self.xAchse = pl.arange(0, 100, 1)
        self.yAchse = pl.array([0]*100)
        self.fig = pl.figure(1)
        self.ax = self.fig.add_subplot(111)
        self.ax.grid(True)
        self.ax.set_title("Realtime Waveform Plot")
        self.ax.set_xlabel("Time")
        self.ax.set_ylabel("Amplitude")
        self.ax.axis([0, 100, -1.5, 1.5])
        self.line1 = self.ax.plot(self.xAchse, self.yAchse, '-')
        self.manager = pl.get_current_fig_manager()

    def RealtimePloter(self, values):
        CurrentXAxis = pl.arange(len(values)-100, len(values), 1)
        self.line1[0].set_data(CurrentXAxis, pl.array(values[-100:]))
        self.ax.axis([CurrentXAxis.min(), CurrentXAxis.max(), -1.5, 1.5])
        self.manager.canvas.draw()


class RealtimeController(object):
    def __init__(self):
        values = [0.0 for x in range(100)]
        Ta = 0.01
        fcos = 3.5
        Konstant = cos(2*pi*fcos*Ta)
        self.model = RealtimeModel(values, Konstant, 0.0, Konstant, 5.0, 1.0)
        self.view = RealtimeView()

    def run(self):
        timer = self.view.fig.canvas.new_timer(interval=50)
        timer.add_callback(self.view.RealtimePloter, self.model.values)
        timer2 = self.view.fig.canvas.new_timer(interval=50)
        timer2.add_callback(self.model.FMwaveformGenerator, ())
        timer.start()
        timer2.start()
        pl.show()


def main():
    app = RealtimeController()
    app.run()
