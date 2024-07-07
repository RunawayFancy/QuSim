import sys
import numpy as np
from vispy import scene, app ; 
app.use_app(backend_name="pyqt5")
from vispy.geometry.rect import Rect;
import multiprocessing as mp
from qusim.Utils.instr import INSTR


try: mp.set_start_method("spawn");
except : pass ; 

colorIndex = [
  "#FF0000", 
  "#2F7FFF", 
  "#00FF00", 
  "#FFFFFF" 
]; 


def move(r:Rect , disp)->Rect : 
  sz = r.size;
  return Rect( ( r.left  + disp[0]*sz[0], r.bottom + disp[1]*sz[1] ), sz ) ;  

def resize(r:Rect , resz)->Rect : 
  sz = r.size;
  nsz =  ( (1 + resz[0]) * sz[0] , (1 + resz[1])* sz[1] ) ; 
  return Rect( ( r.center[0] - nsz[0]/2, r.center[1] - nsz[1]/2 ), nsz ) ;  


class WaveTrace():
  def __init__(self):
    self.pool = {};  

  def run(self , fifo_ ) :
    # self.m1 = 0.1;
    # self.m2 = 0.3;
    # self.r1 = 0.1;
    # self.r2 = 0.3;
    # self.fifo = fifo_ ; 
    self.canvas = scene.SceneCanvas(keys='interactive', show=True)
    # self.grid = self.canvas.central_widget.add_grid(margin=10) ; 
    # self.grid.spacing = 0
    # self.gridLines = scene.visuals.GridLines(color=(1,1,1,1));
    # self.title = scene.Label("Wave Monitor (F to view all , Shift+C to clear all)", color='w')
    # self.title.height_max = 40 ;
    # self.grid.add_widget(self.title, row=0, col=0, col_span=2)
    # self.yaxis = scene.AxisWidget(orientation='left',
    #                          axis_label='Amplitudes',
    #                          axis_font_size=12,
    #                          axis_label_margin=30,
    #                          tick_label_margin=5)
    
    # self.yaxis.width_max = 50
    # self.grid.add_widget(self.yaxis, row=1, col=0)
    # self.xaxis = scene.AxisWidget(orientation='bottom',
    #                          axis_label='Time',
    #                          axis_font_size=12,
    #                          axis_label_margin=20,
    #                          tick_label_margin=10)
    
    # self.xaxis.height_max = 50
    # self.grid.add_widget(self.xaxis, row=2, col=1)
    
    # self.right_padding = self.grid.add_widget(row=1, col=2, row_span=1)
    # self.right_padding.width_max = 20
    
    # self.view = self.grid.add_view(row=1, col=1, border_color='white')
    # self.view.add(self.gridLines);
    # self.view.camera = 'panzoom'
    # self.xaxis.link_view(self.view) ;
    # self.yaxis.link_view(self.view) ;
    
    # self.tmr = app.Timer(0.0166666666666 , connect = self.update  ,start = True) ; 
    # self.canvas.connect(self.on_key_press);

    print('1111')
    
    app.run() ;


  def on_key_press(self, event):
    if(event.text == 'F' or event.text == 'f') :
      #print("F is pressed") ; 
      ts = [];
      bs = [] ;
      ls = []
      rs = []
      hasdata = False; 
      for nn,ps in self.pool.items():
        for p in ps[0] : 
          l,r = p.bounds(0);
          b,t = p.bounds(1);
          ls.append(l) ; rs.append(r) ; 
          ts.append(t) ; bs.append(b) ; 
          hasdata = True ; 
     
      if(hasdata) :
        lm = np.min(ls) ; 
        bm = np.min(bs) ; 
        tM = np.max(ts) ; 
        rM = np.max(rs) ; 
        w = rM - lm  ;
        h  = tM - bm ; 
        rat = 0.1;  
        self.view.camera.rect = Rect((lm - rat*w, bm-rat*h  , max((1+2*rat)*w , 0 ),max((1+2*rat)*h , 0 ))) ;  
        self.view.update();
     
    elif(event.text=='a'):
      self.view.camera.rect = move(self.view.camera.rect,(-self.m1,0))  ;
    elif(event.text=='s'):
      self.view.camera.rect = move(self.view.camera.rect,(0,-self.m1))  ;
    elif(event.text=='w'):
      self.view.camera.rect = move(self.view.camera.rect,(0,self.m1))  ;
    elif(event.text=='d'):
      self.view.camera.rect = move(self.view.camera.rect,(self.m1,0))  ;
    elif(event.text=='A'):
      self.view.camera.rect = move(self.view.camera.rect,(-self.m2,0))  ;
    elif(event.text=='S'):
      self.view.camera.rect = move(self.view.camera.rect,(0,-self.m2))  ;
    elif(event.text=='W'):
      self.view.camera.rect = move(self.view.camera.rect,(0,self.m2))  ;
    elif(event.text=='D'):
      self.view.camera.rect = move(self.view.camera.rect,(self.m2,0))  ;
    elif(event.text=='='):
      self.view.camera.rect = resize(self.view.camera.rect,(0,-self.r1))  ;
    elif(event.text=='-'):
      self.view.camera.rect = resize(self.view.camera.rect,(0,self.r1))  ;
    elif(event.text=='+'):
      self.view.camera.rect = resize(self.view.camera.rect,(-self.r1,0))  ;
    elif(event.text=='_'):
      self.view.camera.rect = resize(self.view.camera.rect,(self.r1,0))  ;
    elif(event.text=='e'):
      self.view.camera.rect = resize(self.view.camera.rect,(-self.r1,-self.r1))  ;
    elif(event.text=='q'):
      self.view.camera.rect = resize(self.view.camera.rect,(self.r1,self.r1))  ;
    elif(event.text=='E'):
      self.view.camera.rect = resize(self.view.camera.rect,(-self.r2,-self.r2))  ;
    elif(event.text=='Q'):
      self.view.camera.rect = resize(self.view.camera.rect,(self.r2,self.r2))  ;
    elif(event.text == 'C'):
      self.clear_all();  

    self.view.update();
       

  def plt(self, v):
    nm =  v[0] ; 
    offs = v[1];
    trcs =  v[2] ; 
    if(nm not in self.pool) :
      text = scene.Text(text=nm , parent=self.view.scene  ,color = "#FFFFFF" ,font_size=16 , anchor_x = 'left' ); 
      text.pos = (trcs[0][-1],offs);
      self.pool[nm] = [[],text] ;
       
    lp = len(self.pool[nm][0]);
    for i,tr in enumerate(trcs[1:]):
      if( i < lp) : 
        self.pool[nm][0][i].set_data(np.vstack((trcs[0],tr+offs)).T) ; 
      else :
        self.pool[nm][0].append( scene.Line( np.vstack((trcs[0],tr+offs)).T , parent=self.view.scene,color = colorIndex[i%4]) ); 
    self.pool[nm][1].pos = (trcs[0][-1], offs );
  
  def update(self , a  = None ):
    print('call update func')
    has = False; 
    while(not self.fifo.empty() ):
      w =  self.fifo.get();
      print(w)
      if hasattr(w,"__maj_type__") and  getattr(w,"__maj_type__") == "CLEARALL" :
        self.clear_all(); 
      elif hasattr(w,"__maj_type__") and  getattr(w,"__maj_type__") == "VIEW" :
        self.view.camera.rect = Rect(w.v)  ;
        self.view.update(); 
      elif(isinstance(w, dict)): 
        for v in w.values() : self.plt(v)
      else:
        self.plt(w);
      has = True;
    if(has) : self.view.update(); 

  def clear_all(self, a = None):
    for nm in self.pool:
      self.pool[nm][1].parent=None;
      for pls in self.pool[nm][0]:
        pls.parent = None ;  
    self.pool.clear();


