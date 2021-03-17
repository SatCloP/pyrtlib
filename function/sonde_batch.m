function batch = sonde_batch(serial_number,flag)
%
% function batch = sonde_batch(serial_number)
%
% interpret Vaisala sonde serial numbers
%
% taken from ARM Sonde site:
%
%                Prior to October, 1994 the Vaisala RS-80 serial
%                number code was:
%
%                DDMMYTTPP, in which
%
%                   DD = day of the month (1-31)
%                   MM = month (1-12) + facility identifier (00, 20,
%                   40, or 80)
%                    Y = last digit of the year     
%                   TT = calibration tray identifier
%                   PP = position in calibration tray
%
%                 More recent radiosonde serial numbers are coded
%
%                 YWWDTTTNN, in which
%
%                    Y = last digit of the year
%                   WW = week number (1-52)
%                    D = day of the week (1-7) Monday=1
%                  TTT = calibration tray identifer
%                   NN = position in calibration tray
%
%                 RS-90 radiosondes (yet to be used operationally by
%                 ARM) are coded
%
%                 YWWDSSSS, in which
%
%                    Y = alphabetic code for the year (T=1998, U=1999,
%                    etc.)
%                   WW = week number (1-52)
%                    D = day of the week (1-7) Monday=1
%                 SSSS = sequence number
%
% DCT 9-11-2000
%