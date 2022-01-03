classdef cougarConstants
   
    properties (Constant=true)
        radian2arcsec = 180*3600/pi;
        radian2arcmin = 180*60/pi;
        arcsec2radian = pi/180/3600;
        arcmin2radian = pi/180/60;
        plank         = 6.62606896e-34; % J/s=m^2kg/s
        c             = 299792458; % light of speed m/s
    end
    
end