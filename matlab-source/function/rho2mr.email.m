Date: Fri, 15 Sep 2000 17:07:16 +0000 (UTC)
From: Dave Turner <dturner@grapple.pnl.gov>
Reply-To: Dave Turner <dave.turner@pnl.gov>
To: Dave Tobin <DaveT@ssec.wisc.edu>
Subject: Conversion Eqs


Hi Dave,

  The conversion equations for you:

   ; rho in g/m3
   ; t   in C
   ; p   in mb
   function rho2w, rho, t, p
     w = rho * (t + 273.16) / (p * 0.3477)
   return,w

   ; w in g/kg
   ; t in C
   ; p in mb
   function w2rho, w, t, p
     rho = w * p * 0.3477 / (t + 273.16)
   return,rho

This was provide by Holger Linne' from Max Planck Institute.

  Have a good trip!
Dave

--
Dave Turner
Pacific Northwest National Laboratory
Currently at the University of Wisconsin-Madison
dave.turner@pnl.gov

