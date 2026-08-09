from qusim.data_view.tracer import *


def run_one(fifo_):
    w = WaveTrace(); 
    w.run(fifo_) ; 
    

class WaveTraceDevice(Instr) :
  def __init__(self) :
    Instr.__init__(self) ;   
    self.q = mp.Queue(1);
    self.put(name='aa', offs=1, traces = np.linspace(0,10,10))
    self.P = mp.Process( target= run_one, args=(self.q , ) ); 
    self.P.run() ;

  
  def put(self , name , offs, traces):
    self.q.put( [name , offs, traces] );

  def clear(self, k = None) :
    self.q.put( trc.CLEARALL() ); 
    
  def __del__(self):
    self.P.kill(); 
  
  W={
    "put":put,  
    "clear":clear
  }