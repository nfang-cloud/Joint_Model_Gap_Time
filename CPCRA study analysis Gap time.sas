
*Import the dataset from lib ddiddc;
data ddiddc;
set ddiddc.ddiddc;
trt=randgrp-1;
gender=gender-1;
hemobl=hemobl-12;
cd4bl=cd4bl/100;
id=seq;

run;


data one;
set ddiddc;
array event_all{20} MAC PCP1 PCP2 PCP3 CANE1 CANE2 CANE3 CMV WAST 
KSV ADC CRYC TB LYMP PML TOXO CRYS OMYC HIST HZ1;
array time_all{20} T2MAC T2PCP1 T2PCP2 T2PCP3 T2CANE1 T2CANE2 T2CANE3 T2CMV T2WAST 
T2KSV T2ADC T2CRYC T2TB T2LYMP T2PML T2TOXO T2CRYS T2OMYC T2HIST T2HZ1;
aa=1;
do i=1 to 20;
	event=event_all{i};
	stoptime=time_all{i}/30;
	output;
end;
run;

* The dataset for recurrent event;
data two;
set one;
if event=1;
run;

proc sort data=two;
by id stoptime;
run;


* The dataset for death event;
data three;
set ddiddc;
by seq;
if first.seq;
stoptime=t2death/30;
event=death*2;		* Set event=2 for death;
aa=1;
run;
* Get the quantiles for death time;
proc univariate data=three noprint;
var stoptime; 
output out=quant_d pctlpts=0 10 20 30 40 50 60 70 80 90 100 pctlpre=qd; 
where event=2;
run;
data quant_d;
set quant_d;
qd0=0;
aa=1;
run;
* Merge the recurrent and death event times;
data four;
set two three;
run;

proc sort data=four;
by id stoptime;
run;
* Calculate the number of recurrent events as a mediator for death;

data four2;
set four;
by id;
retain last_stop nevent;
if first.id then do;
	nevent=0;
	start=0;
	stop=stoptime;
	last_stop=stoptime;
end;
else do;
	nevent=nevent+1;
	start=last_stop;
	stop=stoptime;
	last_stop=stoptime;

end;
gap=stoptime-start;
FstOccur=0;
ReOccur=0;

IF nevent=0 THEN FstOccur=1;
IF nevent>=1 THEN ReOccur=1;
run;
* Get the quantiles for recurrent events;
proc univariate data=four2 noprint;
var stoptime; 
output out=quant0_r pctlpts=0 10 20 30 40 50 60 70 80 90 100 pctlpre=qr0_; 
where FstOccur=1 and event=1;
run;
proc univariate data=four2 noprint;
var stoptime; 
output out=quant1_r pctlpts=0 10 20 30 40 50 60 70 80 90 100 pctlpre=qr1_; 
where ReOccur=1 and event=1;
run;

data quant0_r;
set quant0_r;
aa=1;
qr0_0=0;
run;
data quant1_r;
set quant1_r;
aa=1;
qr1_0=0;
run;

data four3;
merge four2 quant0_r quant1_r quant_d;
by aa;
run;

data five;
set four3;
array quant0_r {4} qr0_0 qr0_30 qr0_60 qr0_100;
array quant1_r {4} qr1_0 qr1_30 qr1_60 qr1_100;
array quant_d {11} qd0 qd10 qd20 qd30 qd40 qd50 qd60 qd70 qd80 qd90 qd100;

array dur0_r {3} dur0_r1-dur0_r3;
array dur1_r {3} dur1_r1-dur1_r3;

array dur_d {10} dur_d1-dur_d10;

array event0_r {3} event0_r1-event0_r3;
array event1_r {3} event1_r1-event1_r3;

array event_d {10} event_d1-event_d10;

array median0_r {3} median0_r1-median0_r3;
array median1_r {3} median1_r1-median1_r3;

array median_d {10} median_d1-median_d10;

do i=1 to 3;
	dur0_r{i}=0;
	dur1_r{i}=0;
    
	event0_r{i}=0;
	event1_r{i}=0;
	
	median0_r{i}=0;
	median1_r{i}=0;
	
end;

do i=1 to 10;
	dur_d{i}=0;
	event_d{i}=0;
	median_d{i}=0;
end;

* For recurrent event;
IF FstOccur=1 THEN DO;
	do i=2 to 4;
		if gap<=quant0_r{i} then do;
			dur0_r{i-1}=gap-quant0_r{i-1};
			median0_r{i-1}=quant0_r{i-1}+dur0_r{i-1}/2; /* Get the median of each interval */
			lastmed0_r=median0_r{i-1};
			i=4;
		end;
        
		else do;
			dur0_r{i-1}=quant0_r{i}-quant0_r{i-1};
			median0_r{i-1}=quant0_r{i-1}+dur0_r{i-1}/2; /* Get the median of each interval */
			lastmed0_r=median0_r{i-1};
		end;
    end;
end;
 
IF ReOccur=1 THEN DO;
	do i=2 to 4;
		if gap<=quant1_r{i} then do;
			dur1_r{i-1}=gap-quant1_r{i-1};
			median1_r{i-1}=quant1_r{i-1}+dur1_r{i-1}/2; /* Get the median of each interval */
			lastmed1_r=median1_r{i-1};
			i=4;
		end;
        
		else do;
			dur1_r{i-1}=quant1_r{i}-quant1_r{i-1};
			median1_r{i-1}=quant1_r{i-1}+dur1_r{i-1}/2; /* Get the median of each interval */
			lastmed1_r=median1_r{i-1};
		end;
   end;
end;


* For recurrent event;
if event=1 then do;
	do i=2 to 4;
		if (gap<=quant0_r{i} AND FstOccur=1) then do;
			event0_r{i-1}=1;
			i=4;
		end;
        if (gap<=quant1_r{i} AND ReOccur=1) then do;
			event1_r{i-1}=1;
			i=4;
		end;

		
	end;

end;


else do; /* If death or censored observation */
	do i=2 to 11;
		if stoptime<=quant_d{i} then do;
			event_d{i-1}=(event=2);
			dur_d{i-1}=stoptime-quant_d{i-1};
			median_d{i-1}=quant_d{i-1}+dur_d{i-1}/2; /* Get the median of each interval */
			lastmed_d=median_d{i-1};
			i=11;
		end;
		else do;
			dur_d{i-1}=quant_d{i}-quant_d{i-1};
			median_d{i-1}=quant_d{i-1}+dur_d{i-1}/2; /* Get the median of each interval */
			lastmed_d=median_d{i-1};
		end;
	end;
end;

run;


* For each of mediator "nevent", we need to calculate the duration in each quantile interval of death time;
data five2;
set five;
array quant {11} qd0 qd10 qd20 qd30 qd40 qd50 qd60 qd70 qd80 qd90 qd100;
array dur {10} dur1-dur10;

last_start=start;
do i=1 to 10;
	dur{i}=0;
end;

do i=2 to 11;
	if stop>quant{i} then do;
		if last_start<= quant{i} then do;
			dur{i-1}=quant{i}-last_start;
			last_start=quant{i};
		end;
	end;
	else do;
		dur{i-1}=stop - last_start;
		i=11;
	end;
end;
run;


**************************************************************************************************************
**************************************************************************************************************
************************************************* Model I ****************************************************;

proc nlmixed data=five2 qpoints=5 MAXITER=5000 corr cov;

parms log0_r1=-2.68 log0_r2=-1.41 log0_r3=-0.39
      log1_r1=-1.47 log1_r2=-0.47 log1_r3=-0.08

	  log_h1=-4.8  log_h2=-3.99 log_h3=-4.21 log_h4=-4.11  log_h5=-3.68 log_h6=-3.68
	  log_h7=-3.15  log_h8=-3.11 log_h9=-3.41 log_h10=-2.90

      betaz=-0.06 betax1=-0.6 betax2=0.4 betax3=-0.07 betax4=-0.1 betax5=-0.5
     
      etaz=-0.40 etax1=0.86 etax2=0.20 etax3=-0.1 etax4=-0.35 etax5=-0.77
	  etam=0.27

      delta1=2.15
	  log_varc=-1.66;
        
array dur {10} dur1-dur10;
array baseh {10} log_h1 log_h2 log_h3 log_h4 log_h5 log_h6 log_h7 log_h8 log_h9 log_h10;

base_haz_r0=exp(log0_r1) * event0_r1 + exp(log0_r2) * event0_r2 + exp(log0_r3) * event0_r3;

base_haz_r1=exp(log1_r1) * event1_r1 + exp(log1_r2) * event1_r2 + exp(log1_r3) * event1_r3;


cum_base_haz_r0=exp(log0_r1) * dur0_r1 + exp(log0_r2) * dur0_r2 + exp(log0_r3) * dur0_r3;

cum_base_haz_r1=exp(log1_r1) * dur1_r1 + exp(log1_r2) * dur1_r2 + exp(log1_r3) * dur1_r3;


mu1= betaz * trt + betax1 * prevoi + betax2 * Gender + 
       betax3 * stratum +  betax4 * hemobl +  betax5 * cd4bl  + vi;			/* for recurrent event first occur */


loglik1=-exp(mu1) *cum_base_haz_r0*FstOccur - exp(mu1)*cum_base_haz_r1*ReOccur;

* Note here the mediator: nevent changes over time, so we need to calculate it for each interval;
* Each record is an interval (start, stop) for nevent;
* dur{k} is the duration of nevent in each of the death quantiles, where nevent is a constant in each dur{k};
sum2=0;
do k=1 to 10;
	/* cumulative baseline hazard for time dependent measure */
sum2=sum2 + exp(baseh{k}) * dur{k} * exp(etam * nevent );
	
end;

mu2= etaz * trt + etax1 * prevoi + etax2 * Gender + 
       etax3 * stratum +  etax4 * hemobl +  etax5 * cd4bl  + delta1 * vi;	/* for death event */
loglik0=-exp(mu2) * sum2;

if event=2 then do; 				/* for death event */
	base_haz_d=exp(log_h1) * event_d1 + exp(log_h2) * event_d2 + exp(log_h3) * event_d3
               + exp(log_h4) * event_d4 + exp(log_h5) * event_d5 + exp(log_h6) * event_d6
               + exp(log_h7) * event_d7 + exp(log_h8) * event_d8 + exp(log_h9) * event_d9
               + exp(log_h10) * event_d10;
	mu4= etaz * trt + etax1 * prevoi + etax2 * Gender + 
       etax3 * stratum +  etax4 * hemobl +  etax5 * cd4bl + etam * nevent + delta1 * vi ;
	
end;

if (event=1 and FstOccur=1) then loglik= log(base_haz_r0)+ mu1 +loglik0 + loglik1;
if (event=1 and ReOccur=1) then loglik= log(base_haz_r1)+ mu1 +loglik0 + loglik1;	/*log likelihood for recurrent event*/

if event=2 then loglik=loglik0 + log(base_haz_d) + mu4 + loglik1;	/*log likelihood for death*/
if event=0 then loglik=loglik0 + loglik1;		/*log likelihood for censoring*/

model id ~ general(loglik);
random vi ~ normal(0,  exp(log_varc)) subject=id;
ods output ParameterEstimates=est1 FitStatistics=fit1 CorrMatParmEst=corr1 CovMatParmEst=cov1;
 
run;



*&out is the path for saving results;
data _null_;
	set cov1;
	file "&out.\covI_realgap.txt";
	put parameter 
		log0_r1 log0_r2 log0_r3
      	log1_r1 log1_r2 log1_r3
      	log_h1  log_h2 log_h3 log_h4  log_h5 log_h6
	  	log_h7  log_h8 log_h9 log_h10
		betaz betax1 betax2 betax3 betax4 betax5
     	etaz etax1 etax2 etax3 etax4 etax5
	  	etam delta1 log_varc;
	run;

data _null_;
	set est1;
	file "&out.\estI_realgap.txt";
	put Parameter Estimate;
	run;

**************************************************************************************************************
**************************************************************************************************************
************************************************* Model II ****************************************************;

proc nlmixed data=five2 qpoints=5 MAXITER=5000 corr cov;

parms log0_r1=-2.68 log0_r2=-1.41 log0_r3=-0.39
      log1_r1=-1.47 log1_r2=-0.47 log1_r3=-0.08

	  log_h1=-4.8  log_h2=-3.99 log_h3=-4.21 log_h4=-4.11  log_h5=-3.68 log_h6=-3.68
	  log_h7=-3.15  log_h8=-3.11 log_h9=-3.41 log_h10=-2.90

      betaz=-0.06 betax1=-0.6 betax2=0.4 betax3=-0.07 betax4=-0.1 betax5=-0.5
     
      etaz=-0.40 etax1=0.86 etax2=0.20 etax3=-0.1 etax4=-0.35 etax5=-0.77
	  etam=0.27

      delta1=2 delta2=0
	  log_varc=-1.66;
        
array dur {10} dur1-dur10;
array baseh {10} log_h1 log_h2 log_h3 log_h4 log_h5 log_h6 log_h7 log_h8 log_h9 log_h10;

base_haz_r0=exp(log0_r1) * event0_r1 + exp(log0_r2) * event0_r2 + exp(log0_r3) * event0_r3;

base_haz_r1=exp(log1_r1) * event1_r1 + exp(log1_r2) * event1_r2 + exp(log1_r3) * event1_r3;


cum_base_haz_r0=exp(log0_r1) * dur0_r1 + exp(log0_r2) * dur0_r2 + exp(log0_r3) * dur0_r3;

cum_base_haz_r1=exp(log1_r1) * dur1_r1 + exp(log1_r2) * dur1_r2 + exp(log1_r3) * dur1_r3;


mu1= betaz * trt + betax1 * prevoi + betax2 * Gender + 
       betax3 * stratum +  betax4 * hemobl +  betax5 * cd4bl  + vi;			/* for recurrent event first occur */


loglik1=-exp(mu1) *cum_base_haz_r0*FstOccur - exp(mu1)*cum_base_haz_r1*ReOccur;

* Note here the mediator: nevent changes over time, so we need to calculate it for each interval;
* Each record is an interval (start, stop) for nevent;
* dur{k} is the duration of nevent in each of the death quantiles, where nevent is a constant in each dur{k};
sum2=0;
do k=1 to 10;
	/* cumulative baseline hazard for time dependent measure */
sum2=sum2 + exp(baseh{k}) * dur{k} * exp(etam * nevent + delta2 * nevent * vi);
	
end;

mu2= etaz * trt + etax1 * prevoi + etax2 * Gender + 
       etax3 * stratum +  etax4 * hemobl +  etax5 * cd4bl  + delta1 * vi;	/* for death event */
loglik0=-exp(mu2) * sum2;

if event=2 then do; 				/* for death event */
	base_haz_d=exp(log_h1) * event_d1 + exp(log_h2) * event_d2 + exp(log_h3) * event_d3
               + exp(log_h4) * event_d4 + exp(log_h5) * event_d5 + exp(log_h6) * event_d6
               + exp(log_h7) * event_d7 + exp(log_h8) * event_d8 + exp(log_h9) * event_d9
               + exp(log_h10) * event_d10;
	mu4= etaz * trt + etax1 * prevoi + etax2 * Gender + 
       etax3 * stratum +  etax4 * hemobl +  etax5 * cd4bl + etam * nevent + delta1 * vi + delta2 * nevent * vi;
	
end;

if (event=1 and FstOccur=1) then loglik= log(base_haz_r0)+ mu1 +loglik0 + loglik1;
if (event=1 and ReOccur=1) then loglik= log(base_haz_r1)+ mu1 +loglik0 + loglik1;	/*log likelihood for recurrent event*/

if event=2 then loglik=loglik0 + log(base_haz_d) + mu4 + loglik1;	/*log likelihood for death*/
if event=0 then loglik=loglik0 + loglik1;		/*log likelihood for censoring*/

model id ~ general(loglik);
random vi ~ normal(0,  exp(log_varc)) subject=id;
ods output ParameterEstimates=est2 FitStatistics=fit2 CorrMatParmEst=corr2 CovMatParmEst=cov2;
 
run;

*&out is the path for saving results;
data _null_;
	set cov2;
	file "&out.\covII_realgap.txt";
	put parameter 
		log0_r1 log0_r2 log0_r3
      	log1_r1 log1_r2 log1_r3
      	log_h1  log_h2 log_h3 log_h4  log_h5 log_h6
	  	log_h7  log_h8 log_h9 log_h10
		betaz betax1 betax2 betax3 betax4 betax5
     	etaz etax1 etax2 etax3 etax4 etax5
	  	etam delta1 delta2 log_varc;
	run;

data _null_;
	set est2;
	file "&out.\estII_realgap.txt";
	put Parameter Estimate;
	run;


