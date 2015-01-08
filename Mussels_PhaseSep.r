# A mussel pattern formation model that follows the Cahn-Hilliard model
# for phase separation. See Liu et al, PNAS 2013, 110 :11905-11910

# First setup of the model
#remove(list=ls()) # Remove all variables from memory

require("fields") # An R package to load

D    = 1;         # 1    - The diffusivity parameter
Beta = 1.89;      # 1.89 - The non-parameterical density-dependence parameter
k1   = 0.1;       # 0.1  - The non-local diffusion term
  
# Simulation parameters
Length=32             #  128 - The length of the simulated landscape, in meters
m=64                  #  128 - # gridcellen
dT=0.01               # 0.25 - timestep
EndTime=1000          #      - Time at which the simulation ends
Numframes=300         #      - Number of frames displayed during the entire simulation

WinWidth=7             #      - Width of the simulation window 
WinHeight=5            #      - Height of the simulation window

NX=NY=m                # Setting the grid dimensions
dx=Length/NX           # The size of a grid cell in the x dimension
dy=Length/NY           # The size of a grid cell in the y dimension

# ------ Function definitions
Laplacian <- function (C, dx, dy) {
  
  NX=dim(C)[2]
  NY=dim(C)[1]
  
  dC=matrix(nrow=NY, ncol=NX, data=0)
  
  dC[2:(NY-1),2:(NX-1)] = 
    (C[2:(NY-1),1:(NX-2)] + C[2:(NY-1),3:NX] - 2*C[2:(NY-1),2:(NX-1)])/dx/dx +
    (C[1:(NY-2),2:(NX-1)] + C[3:NY,2:(NX-1)] - 2*C[2:(NY-1),2:(NX-1)])/dy/dy
  
  return(dC)
  
}

GradientX <- function (C, dx, dy) {
  
  NX=dim(C)[2]
  NY=dim(C)[1]
  
  dC=matrix(nrow=NY, ncol=NX, data=0)
  
  dC[2:(NY-1),2:(NX-1)] = 
    - (C[2:(NY-1),1:(NX-2)] - C[2:(NY-1),3:NX] )/dx/2

  return(dC)
  
}

GradientY <- function (C, dx, dy) {
  
  NX=dim(C)[2]
  NY=dim(C)[1]
  
  dC=matrix(nrow=NY, ncol=NX, data=0)
  
  dC[2:(NY-1),2:(NX-1)] = 
    - (C[1:(NY-2),2:(NX-1)] - C[3:NY,2:(NX-1)] )/dy/2
  
  return(dC)
  
}

# Initialisation: declaring variable and setting up initial values
# All arrays of dimension m x m

C = matrix(nrow=NY,ncol=NX, data=runif(NY*NX)*1.3)
C = matrix(nrow=NY,ncol=NX, data=0)
X = matrix(nrow=15,ncol=15, data=(1:15)*4) 
Y = matrix(nrow=15,ncol=15, data=(1:15)*4, byrow=TRUE) 
C[X,Y]=1.3;
dC = matrix(nrow=NY,ncol=NX, data=0)
J2 = J3 = J4 = matrix(nrow=NY-2,ncol=NX-2, data=0)

RecordTimes=(0:Numframes)*EndTime/Numframes
Step=EndTime^(1/(Numframes-1))
RecordTimes=Step^(0:(Numframes-1))

AllData=array(data=0,dim=c(m,m,length(RecordTimes)))

Time  =  1        ; # Begin time 
jj    =  1        ; # The counter needed for recording data during the run

# ------- Setting up the figure ------------------------------------------
  
## Open a graphics window (Darwin stands for a Mac computer)
if (Sys.info()["sysname"]=="Darwin"){
    quartz(width=WinWidth, height=WinHeight, title="Mussel pattern formation")
  } else
    win.graph(width = WinWidth, height = WinHeight)

par(mar=c(3, 4, 2, 6) + 0.1)

# ------------ The simulation loop ---------------------------------------

print(system.time(
while (Time<=EndTime){  # Here the time loop starts 

  J0 = (C^2 - Beta*C + 1)*(3*C^2 - 2*Beta*C + 1 )
  
  J1X = J0 * GradientX(C,dx,dy)
  J2X = D*GradientX(J1X,dx,dy)
  
  J1Y = J0 * GradientY(C,dx,dy)
  J2Y = D*GradientY(J1Y,dx,dy)
  
  J3 = Laplacian(C,dx,dy)
  J4 = D * k1 * Laplacian(J3,dx,dy)
  
  dC = J2X + J2Y - J4
  
  
    # Graphic representation of the model every now and then
    if ((Time+dT/100)>=RecordTimes[jj]){   # The dT/100 is to avoid inaccuracy problems
    
      image.plot(C, zlim=c(0,1.4), xaxt="n", yaxt="n", useRaster=TRUE,
                 col = topo.colors(255),asp=1, bty="n",
                 legend.shrink = 0.99, legend.width = 2)  
      title("Mussel density")   
    
      mtext(text=paste("Time : ",sprintf("%1.0f",Time),
                       "of" ,sprintf("%1.0f",EndTime), "Minutes"), 
                       side=1, adj=0.5, line=1.5, cex=1)
    
      dev.flush()
      dev.hold()
    
      AllData[,,jj]=C
    
      jj=jj+1 # Increasing the Recorder counter
    
    } 
      
    C = C + dC*dT;
    
    # Circular boundary conditions
    
    C[,1:2]=C[,(NX-3):(NX-2)]
    C[1:2,]=C[(NY-3):(NY-2),]
    C[,(NX-1):NX]=C[,3:4]
    C[(NY-1):NY,]=C[3:4,]
    
    Time=Time+dT;  # Incrementing time by dT
  
    #pause
    
  }

))

save(file='Mussel_PDE.rda',AllData,RecordTimes)

system('say Finished')
