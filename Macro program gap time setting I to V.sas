proc datasets lib=work kill;run;quit;
DM 'log; clear;out;clear;odsresults; clear;';

***********************************************************************************
*************************Setting I***************************************************;

%macro jointmodel1(Data=, parms=, mod=);

%do ii=1 %to 200;
title "Iteration &ii.";
*Import the dataset;
*&InD is the path for import the datasets;
PROC IMPORT OUT= WORK.sim2 
            DATAFILE= "&InD.\&Data.&ii..csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;
*Import the jump points;
data quant0_r;
infile "&InD.\quant_r.txt";
input qr0_min qr0_1 qr0_max aa;
run;

data quant1_r;
infile "&InD.\quant_r.txt";
input qr1_min qr1_1 qr1_max aa;
run;

data quant2_r;
infile "&InD.\quant_r.txt";
input qr2_min qr2_1 qr2_max aa;
run;

data quant_d;
infile "&InD.\quant_d.txt";
input qd_min qd_1 qd_2 qd_3 qd_max aa;
run;

data four;
set sim2;
aa=1;
run;

proc sort data=four;
by id stoptime;
run;

* Calculate the number of recurrent events as a mediator for death and gap time;
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
ThrOccur=0;
IF nevent=0 THEN FstOccur=1;
IF nevent=1 THEN ReOccur=1;
IF nevent>=2 THEN ThrOccur=1;
run;

data four3;
merge four2 quant_d quant0_r quant1_r quant2_r;
by aa;
run;

* Calculate the duration in each quantile interval, together with the indicator of event in each interval;
data five;
set four3;
array quant0_r {3} qr0_min qr0_1 qr0_max;
array quant1_r {3} qr1_min qr1_1 qr1_max;
array quant2_r {3} qr2_min qr2_1 qr2_max;
array quant_d {4} qd_min qd_1 qd_2 qd_max;

array dur0_r {2} dur0_r1-dur0_r2;
array dur1_r {2} dur1_r1-dur1_r2;
array dur2_r {2} dur2_r1-dur2_r2;
array dur_d {3} dur_d1-dur_d3;

array event0_r {2} event0_r1-event0_r2;
array event1_r {2} event1_r1-event1_r2;
array event2_r {2} event2_r1-event2_r2;
array event_d {3} event_d1-event_d3;

array median0_r {2} median0_r1-median0_r2;
array median1_r {2} median1_r1-median1_r2;
array median2_r {2} median2_r1-median2_r2;
array median_d {3} median_d1-median_d3;

do i=1 to 2;
	dur0_r{i}=0;
	dur1_r{i}=0;
    dur2_r{i}=0;
	event0_r{i}=0;
	event1_r{i}=0;
	event2_r{i}=0;
	median0_r{i}=0;
	median1_r{i}=0;
	median2_r{i}=0;
end;

do i=1 to 3;
	dur_d{i}=0;
	event_d{i}=0;
	median_d{i}=0;
end;

* For recurrent event;
IF FstOccur=1 THEN DO;
	do i=2 to 3;
		if gap<=quant0_r{i} then do;
			dur0_r{i-1}=gap-quant0_r{i-1};
			median0_r{i-1}=quant0_r{i-1}+dur0_r{i-1}/2; /* Get the median of each interval */
			lastmed0_r=median0_r{i-1};
			i=3;
		end;
        
		else do;
			dur0_r{i-1}=quant0_r{i}-quant0_r{i-1};
			median0_r{i-1}=quant0_r{i-1}+dur0_r{i-1}/2; /* Get the median of each interval */
			lastmed0_r=median0_r{i-1};
		end;
    end;
end;
 
IF ReOccur=1 THEN DO;
	do i=2 to 3;
		if gap<=quant1_r{i} then do;
			dur1_r{i-1}=gap-quant1_r{i-1};
			median1_r{i-1}=quant1_r{i-1}+dur1_r{i-1}/2; /* Get the median of each interval */
			lastmed1_r=median1_r{i-1};
			i=3;
		end;
        
		else do;
			dur1_r{i-1}=quant1_r{i}-quant1_r{i-1};
			median1_r{i-1}=quant1_r{i-1}+dur1_r{i-1}/2; /* Get the median of each interval */
			lastmed1_r=median1_r{i-1};
		end;
   end;
end;

IF ThrOccur=1 THEN DO;
	do i=2 to 3;
		if gap<=quant2_r{i} then do;
			dur2_r{i-1}=gap-quant2_r{i-1};
			median2_r{i-1}=quant2_r{i-1}+dur2_r{i-1}/2; /* Get the median of each interval */
			lastmed2_r=median2_r{i-1};
			i=3;
		end;
        
		else do;
			dur2_r{i-1}=quant2_r{i}-quant2_r{i-1};
			median2_r{i-1}=quant2_r{i-1}+dur2_r{i-1}/2; /* Get the median of each interval */
			lastmed2_r=median2_r{i-1};
		end;
   end;
end;
* For recurrent event;
if event=1 then do;
	do i=2 to 3;
		if (gap<=quant0_r{i} AND FstOccur=1) then do;
			event0_r{i-1}=1;
			i=3;
		end;
        if (gap<=quant1_r{i} AND ReOccur=1) then do;
			event1_r{i-1}=1;
			i=3;
		end;

		if (gap<=quant2_r{i} AND ThrOccur=1) then do;
			event2_r{i-1}=1;
			i=3;
		end;
	end;

end;


else do; /* If death or censored observation */
	do i=2 to 4;
		if stoptime<=quant_d{i} then do;
			event_d{i-1}=(event=2);
			dur_d{i-1}=stoptime-quant_d{i-1};
			median_d{i-1}=quant_d{i-1}+dur_d{i-1}/2; /* Get the median of each interval */
			lastmed_d=median_d{i-1};
			i=4;
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
array quant {4} qd_min qd_1 qd_2 qd_max;
array dur {3} dur1-dur3;

last_start=start;
do i=1 to 3;
	dur{i}=0;
end;

do i=2 to 4;
	if stop>quant{i} then do;
		if last_start<= quant{i} then do;
			dur{i-1}=quant{i}-last_start;
			last_start=quant{i};
		end;
	end;
	else do;
		dur{i-1}=stop - last_start;
		i=4;
	end;
end;
run;	

proc nlmixed data=five2 qpoints=5 MAXITER=5000 corr cov;

parms &parms.;
        
array dur {3} dur1-dur3;
array baseh {3} log_h1 log_h2 log_h3;

base_haz_r0=exp(log0_r1) * event0_r1 + exp(log0_r2) * event0_r2;

base_haz_r1=exp(log1_r1) * event1_r1 + exp(log1_r2) * event1_r2;

base_haz_r2=exp(log2_r1) * event2_r1 + exp(log2_r2) * event2_r2;

cum_base_haz_r0=exp(log0_r1) * dur0_r1 + exp(log0_r2) * dur0_r2;

cum_base_haz_r1=exp(log1_r1) * dur1_r1 + exp(log1_r2) * dur1_r2;

cum_base_haz_r2=exp(log2_r1) * dur2_r1 + exp(log2_r2) * dur2_r2;

mu1_1= betaz1 * X1 + betax1 * X2 + vi;			/* for recurrent event first occur */
mu1_2= betaz2 * X1 + betax2 * X2 + vi;			/* for recurrent event reoccur */
mu1_3= betaz3 * X1 + betax3 * X2 + vi;			/* for recurrent event reoccur */

loglik1=-exp(mu1_1) *cum_base_haz_r0*FstOccur - exp(mu1_2)*cum_base_haz_r1*ReOccur- exp(mu1_3)*cum_base_haz_r2*ThrOccur;

sum2=0;
do k=1 to 3;
/* cumulative baseline hazard for time dependent measure */
sum2=sum2 + exp(baseh{k}) * dur{k} * exp(etam * nevent );
	
end;

mu2= etaz * X1 + etax * X2 + delta1 * vi;	/* for death event */
loglik0=-exp(mu2) * sum2;

if event=2 then do; 				/* for death event */
	base_haz_d=exp(log_h1) * event_d1 + exp(log_h2) * event_d2 + exp(log_h3) * event_d3;
	mu4= etaz * X1 + etax * X2 + etam * nevent + delta1 * vi ;
	
end;

if (event=1 and FstOccur=1) then loglik= log(base_haz_r0)+ mu1_1 +loglik0 + loglik1; /*log likelihood for first recurrent event*/
if (event=1 and ReOccur=1) then loglik= log(base_haz_r1)+ mu1_2 +loglik0 + loglik1;	/*log likelihood for recuring recurrent event*/
if (event=1 and ThrOccur=1) then loglik= log(base_haz_r2)+ mu1_3 +loglik0 + loglik1;/*log likelihood for 3rd recurrent event*/

if event=2 then loglik=loglik0 + log(base_haz_d) + mu4 + loglik1;	/*log likelihood for death*/
if event=0 then loglik=loglik0 + loglik1;		/*log likelihood for censoring*/

model id ~ general(loglik);
random vi ~ normal(0,  exp(log_varc)) subject=id;
ods output ParameterEstimates=est&mod.&ii. FitStatistics=fit&mod.&ii. CorrMatParmEst=corr&mod.&ii. CovMatParmEst=cov&mod.&ii.;
 
run;
%end;

%mend;
%Let parms_list1=log0_r1=1.1 log0_r2=1.16 
      log1_r1=1.28 log1_r2=1.50 
      log2_r1=1.44 log2_r2=1.61 
	  log_h1=-0.1  log_h2=0 log_h3=0.1
	  betax1=0.5 betax2=0.5 betax3=0.5 
      betaz1=-0.5 betaz2=-0.4 betaz3=-0.3
	  etaz=-1 etax=1
	  etam=0.25
      delta1=1
	  log_varc=0;

%jointmodel1(data=sim1_gapdata,parms=&parms_list1., mod=1);

/* &out is the saving file path*/
*Output the covariance;
%macro outcov(mod=,filename=,delta2=);
%do ii=1 %to 200;
data _null_;
	set cov&mod&ii;
	file "&out.\&filename&ii..txt";
	put parameter log0_r1 log0_r2 
       log1_r1 log1_r2 
       log2_r1 log2_r2
	   log_h1 log_h2  log_h3
		betax1 betax2 betax3
        betaz1 betaz2 betaz3
		etaz etax
		etam
        delta1 &delta2
		log_varc
	   ;
	run;
%end;
%mend;

*Output the estimation;
%macro outest(mod=,filename=);
%do ii=1 %to 200;
data _null_;
	set est&mod&ii;
	file "&out.\&filename&ii..txt";
	put Parameter Estimate;
	run;
%end;
%mend;

%outest(mod=1,filename=covIgap,delta2=);
%outcov(mod=1,filename=estIgap);


***********************************************************************************
*************************Setting II***************************************************;


%macro jointmodel2(Data=, parms=, mod=);

%do ii=1 %to 200;
title "Iteration &ii.";
*&InD is the path for import the datasets;
PROC IMPORT OUT= WORK.sim2 
            DATAFILE= "&InD.\&Data.&ii..csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

data quant0_r;
infile "&InD.\quant_r.txt";
input qr0_min qr0_1 qr0_max aa;
run;

data quant1_r;
infile "&InD.\quant_r.txt";
input qr1_min qr1_1 qr1_max aa;
run;

data quant2_r;
infile "&InD.\quant_r.txt";
input qr2_min qr2_1 qr2_max aa;
run;

data quant_d;
infile "&InD.1\quant_d.txt";
input qd_min qd_1 qd_2 qd_3 qd_max aa;
run;

data four;
set sim2;
aa=1;
run;

proc sort data=four;
by id stoptime;
run;

* Calculate the number of recurrent events as a mediator for death and gap time;
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
ThrOccur=0;
IF nevent=0 THEN FstOccur=1;
IF nevent=1 THEN ReOccur=1;
IF nevent>=2 THEN ThrOccur=1;
run;

data four3;
merge four2 quant_d quant0_r quant1_r quant2_r;
by aa;
run;

* Calculate the duration in each quantile interval, together with the indicator of event in each interval;
data five;
set four3;
array quant0_r {3} qr0_min qr0_1 qr0_max;
array quant1_r {3} qr1_min qr1_1 qr1_max;
array quant2_r {3} qr2_min qr2_1 qr2_max;
array quant_d {4} qd_min qd_1 qd_2 qd_max;

array dur0_r {2} dur0_r1-dur0_r2;
array dur1_r {2} dur1_r1-dur1_r2;
array dur2_r {2} dur2_r1-dur2_r2;
array dur_d {3} dur_d1-dur_d3;

array event0_r {2} event0_r1-event0_r2;
array event1_r {2} event1_r1-event1_r2;
array event2_r {2} event2_r1-event2_r2;
array event_d {3} event_d1-event_d3;

array median0_r {2} median0_r1-median0_r2;
array median1_r {2} median1_r1-median1_r2;
array median2_r {2} median2_r1-median2_r2;
array median_d {3} median_d1-median_d3;

do i=1 to 2;
	dur0_r{i}=0;
	dur1_r{i}=0;
    dur2_r{i}=0;
	event0_r{i}=0;
	event1_r{i}=0;
	event2_r{i}=0;
	median0_r{i}=0;
	median1_r{i}=0;
	median2_r{i}=0;
end;

do i=1 to 3;
	dur_d{i}=0;
	event_d{i}=0;
	median_d{i}=0;
end;

* For recurrent event;
IF FstOccur=1 THEN DO;
	do i=2 to 3;
		if gap<=quant0_r{i} then do;
			dur0_r{i-1}=gap-quant0_r{i-1};
			median0_r{i-1}=quant0_r{i-1}+dur0_r{i-1}/2; /* Get the median of each interval */
			lastmed0_r=median0_r{i-1};
			i=3;
		end;
        
		else do;
			dur0_r{i-1}=quant0_r{i}-quant0_r{i-1};
			median0_r{i-1}=quant0_r{i-1}+dur0_r{i-1}/2; /* Get the median of each interval */
			lastmed0_r=median0_r{i-1};
		end;
    end;
end;
 
IF ReOccur=1 THEN DO;
	do i=2 to 3;
		if gap<=quant1_r{i} then do;
			dur1_r{i-1}=gap-quant1_r{i-1};
			median1_r{i-1}=quant1_r{i-1}+dur1_r{i-1}/2; /* Get the median of each interval */
			lastmed1_r=median1_r{i-1};
			i=3;
		end;
        
		else do;
			dur1_r{i-1}=quant1_r{i}-quant1_r{i-1};
			median1_r{i-1}=quant1_r{i-1}+dur1_r{i-1}/2; /* Get the median of each interval */
			lastmed1_r=median1_r{i-1};
		end;
   end;
end;

IF ThrOccur=1 THEN DO;
	do i=2 to 3;
		if gap<=quant2_r{i} then do;
			dur2_r{i-1}=gap-quant2_r{i-1};
			median2_r{i-1}=quant2_r{i-1}+dur2_r{i-1}/2; /* Get the median of each interval */
			lastmed2_r=median2_r{i-1};
			i=3;
		end;
        
		else do;
			dur2_r{i-1}=quant2_r{i}-quant2_r{i-1};
			median2_r{i-1}=quant2_r{i-1}+dur2_r{i-1}/2; /* Get the median of each interval */
			lastmed2_r=median2_r{i-1};
		end;
   end;
end;
* For recurrent event;
if event=1 then do;
	do i=2 to 3;
		if (gap<=quant0_r{i} AND FstOccur=1) then do;
			event0_r{i-1}=1;
			i=3;
		end;
        if (gap<=quant1_r{i} AND ReOccur=1) then do;
			event1_r{i-1}=1;
			i=3;
		end;

		if (gap<=quant2_r{i} AND ThrOccur=1) then do;
			event2_r{i-1}=1;
			i=3;
		end;
	end;

end;


else do; /* If death or censored observation */
	do i=2 to 4;
		if stoptime<=quant_d{i} then do;
			event_d{i-1}=(event=2);
			dur_d{i-1}=stoptime-quant_d{i-1};
			median_d{i-1}=quant_d{i-1}+dur_d{i-1}/2; /* Get the median of each interval */
			lastmed_d=median_d{i-1};
			i=4;
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
array quant {4} qd_min qd_1 qd_2 qd_max;
array dur {3} dur1-dur3;

last_start=start;
do i=1 to 3;
	dur{i}=0;
end;

do i=2 to 4;
	if stop>quant{i} then do;
		if last_start<= quant{i} then do;
			dur{i-1}=quant{i}-last_start;
			last_start=quant{i};
		end;
	end;
	else do;
		dur{i-1}=stop - last_start;
		i=4;
	end;
end;
run;	

* Model II in Biometrics paper;
proc nlmixed data=five2 qpoints=5 MAXITER=5000 corr cov;
parms &parms. ;
* dur{k} is the duration of nevent in each of the death quantiles, where nevent is a constant in each dur{k};        
array dur {3} dur1-dur3;
array baseh {3} log_h1 log_h2 log_h3;

base_haz_r0=exp(log0_r1) * event0_r1 + exp(log0_r2) * event0_r2;

base_haz_r1=exp(log1_r1) * event1_r1 + exp(log1_r2) * event1_r2;

base_haz_r2=exp(log2_r1) * event2_r1 + exp(log2_r2) * event2_r2;

cum_base_haz_r0=exp(log0_r1) * dur0_r1 + exp(log0_r2) * dur0_r2;

cum_base_haz_r1=exp(log1_r1) * dur1_r1 + exp(log1_r2) * dur1_r2;

cum_base_haz_r2=exp(log2_r1) * dur2_r1 + exp(log2_r2) * dur2_r2;

mu1_1= betaz1 * X1 + betax1 * X2 + vi;			/* for recurrent event first occur */
mu1_2= betaz2 * X1 + betax2 * X2 + vi;			/* for recurrent event reoccur */
mu1_3= betaz3 * X1 + betax3 * X2 + vi;			/* for recurrent event reoccur */

loglik1=-exp(mu1_1) *cum_base_haz_r0*FstOccur - exp(mu1_2)*cum_base_haz_r1*ReOccur- exp(mu1_3)*cum_base_haz_r2*ThrOccur;

sum2=0;
do k=1 to 3;
	/* cumulative baseline hazard for time dependent measure */
sum2=sum2 + exp(baseh{k}) * dur{k} * exp(etam * nevent + delta2 * nevent * vi);
	
end;

mu2= etaz * X1 + etax * X2 + delta1 * vi;	/* for death event */
loglik0=-exp(mu2) * sum2;

if event=2 then do; 				/* for death event */
	base_haz_d=exp(log_h1) * event_d1 + exp(log_h2) * event_d2 + exp(log_h3) * event_d3;
	mu4= etaz * X1 + etax * X2 + etam * nevent + delta1 * vi + delta2 * nevent * vi;
	
end;

if (event=1 and FstOccur=1) then loglik= log(base_haz_r0)+ mu1_1 +loglik0 + loglik1;
if (event=1 and ReOccur=1) then loglik= log(base_haz_r1)+ mu1_2 +loglik0 + loglik1;	/*log likelihood for recurrent event*/
if (event=1 and ThrOccur=1) then loglik= log(base_haz_r2)+ mu1_3 +loglik0 + loglik1;	/*log likelihood for recurrent event*/

if event=2 then loglik=loglik0 + log(base_haz_d) + mu4 + loglik1;	/*log likelihood for death*/
if event=0 then loglik=loglik0 + loglik1;		/*log likelihood for censoring*/

model id ~ general(loglik);
random vi ~ normal(0,  exp(log_varc)) subject=id;
ods output ParameterEstimates=est&mod.&ii. FitStatistics=fit&mod.&ii. CorrMatParmEst=corr&mod.&ii. CovMatParmEst=cov&mod.&ii.;
 
run;
%end;

%mend;
%Let parms_list2=log0_r1=1.1 log0_r2=1.16 
      log1_r1=1.28 log1_r2=1.50 
      log2_r1=1.44 log2_r2=1.61 
	  log_h1=-0.1  log_h2=0 log_h3=0.1
	  betax1=0.5 betax2=0.5 betax3=0.5 
      betaz1=-0.5 betaz2=-0.4 betaz3=-0.3
	  etaz=-1 etax=1
	  etam=0.25
      delta1=1
	  delta2=0.15
	  log_varc=0;
%jointmodel2(Data=sim2_gapdata, parms=&parms_list2.,mod=2);
*Export the estimation and covariance results;
%outcov(mod=2,filename=covIIgap,delta2=delta2);
%outest(mod=2,filename=estIIgap);

***********************************************************************************
*************************Setting III***************************************************;

*Datasets sim3_alldata generating using log gamma distribution;
%jointmodel2(Data=sim3_gapdata,parms=&parms_list2.,mod=3);
*Export the estimation and covariance results;
%outcov(mod=3,filename=covIIIgap,delta2=delta2);
%outest(mod=3,filename=estIIIgap);


**************************************************************************************
*************************Setting IV***************************************************;

*Fit model without the interaction term while the data are generated with the interaction term;
%jointmodel1(Data=sim2_gapdata,parms=&parms_list1.,mod=4);
*Export the estimation and covariance results;
%outcov(mod=4,filename=covIVgap,delta2=);
%outest(mod=4,filename=estIVgap);

**************************************************************************************
*************************Setting V***************************************************;
%Let parms_list5= log0_r1=0 log0_r2=1.1 log1_r1=1.6 log2_r1=1.8 
	  log_h1=-0.1  log_h2=0 log_h3=0.1
	  betax1=0.5 betax2=0.5 betax3=0.5 
      betaz1=-0.5 betaz2=-0.4 betaz3=-0.3
	  etaz=-1 etax=1
	  etam=0.25
      delta1=1
	  log_varc=0;

%jointmodel1(Data=sim5_gapdata, parms=&parms_list5.,mod=5);
%outcov(mod=5,filename=covVgap,delta2=);
%outest(mod=5,filename=estVgap);
