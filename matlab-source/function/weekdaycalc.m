%This program will calculate the day for any date given by the user.
%The calender works from the year 1582 to the year infinity.
%It is based on the fact that weekdays repeat itself every 400 years.
%A special feature of this calender is that it automatically takes into
%account whether the year is a leapyear or not. 
%
%Just type weekdaycalc to start

fprintf('\n Welcome to the infinite Calender.\n')


theyear=input('\n Please enter the year(example.1999,1932,2056,etc.):    ');


while isempty(theyear)|theyear<1583
   errordlg('Please enter a year(The year should be after 1583 A.D. for reliable results).')
   theyear=input('\nPlease enter the year(example.1999,1932,2056,etc.):    ');



end
   themonth=input('\n Please enter the month:          ');

while isempty(themonth)| themonth>=13

    errordlg('Please enter a number between 1 and 12 to specify the month.')
themonth=input('\n Please enter the month:          ');

end


  theday=input('\n Please enter the date:      ');
  while  isempty(theday)|theday>31|theday<1
     errordlg('Please enter the day between 1 and 31.')
     theday=input('\n Please enter the day:      ');
         end
leapcheck=rem(theyear,400);
R=rem((rem((rem(theyear,1000)),100)),28);
if R==0|R==6|R==17|R==23 
    year=0;
    elseif R==1|R==7|R==12|R==18 
       year=1;
    elseif R==2|R==13|R==19|R==24
       year=2;
elseif R==3|R==8|R==14|R==25
   year=3;
elseif   R==9|R==15|R==20|R==26
     year=4;
elseif R==4|R==10|R==21|R==27
    year=5;
elseif R==5|R==11|R==16|R==22
    year=6;
    else end

C=rem(((theyear-(rem((rem(theyear,1000)),100)))/100),4);
if C==0
   greg=0;
elseif C==1
   greg==5;
elseif C==2
   greg=3;
elseif C==3
   greg=1;
   else end


Day=rem(theday,7);


if themonth==01 & leapcheck==0
   month=6;
elseif themonth==01 & leapcheck~=0
   month=0;
   
elseif themonth==02 & leapcheck==0
month=2;   
elseif themonth==02  & leapcheck~=0
   month=3;
elseif themonth==03|themonth==11
   month=3;
elseif themonth==04|themonth==7
   month=6;
elseif themonth==5
   month=1;
elseif themonth==6
   month=4;
elseif themonth==8
   month=2;
elseif themonth==9|themonth==12
   month=5;
elseif themonth==10
   month=0;
else end

preday=rem((Day+month+year+greg),7);

if preday==0
   fprintf('\n The day of the week is Saturday.\n');
   
   elseif preday==1
   fprintf('\n The day of the week is Sunday.\n');
elseif preday==2
   fprintf('\n The day of the week is Monday.\n');
elseif preday==3
   fprintf('\n The day of the week is Tuesday.\n');
elseif preday==4
   fprintf('\n The day of the week is Wednesday.\n');
elseif preday==5
   fprintf('\n The day of the week is Thursday.\n');
elseif preday==6
   fprintf('\n The day of the week is Friday.\n');
else end

%Made by Manu Mital
