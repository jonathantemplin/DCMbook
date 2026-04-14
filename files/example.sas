* Output options: all can be turned on/off by adding or removing the 'NO':
  page number, date, centering, or page breaks, page length and width;
  OPTIONS nonumber nodate nocenter formdlim=' ' pagesize=MAX linesize=MAX;
* Macro debugging options (add NO in front to turn off);
  OPTIONS nosymbolgen nomerror nomprint nomlogic nospool; 
* Options for calling in and out Mplus: 
  NOXWAIT=return to SAS automatically, NOXSYNC=SAS and Mplus can run at same time;
  OPTIONS NOXWAIT XSYNC;
* Eliminate SAS default titles and names of tables in output (TRACE ON to show);
  TITLE; ODS TRACE OFF;

******************************************************************************
**********       		MACRO VARIABLE DEFINITIONS:			   	      ********
**********		USER NEEDS TO CHANGE THESE VALUES PER ANALYSIS		  ********
******************************************************************************;

* Defining needed macro variables as global;
%GLOBAL macroloc filesave filename saslibname Qname dataname IDname 
		itemstem itemlist numitem ordervar maxitemorder
		attstem attcat numatt numclass structon structorder loosen processors;

* Location of original data files - CHANGE ALL OF THEM;
* Permanent SAS library; 					LIBNAME folder 		"C:\ECPE";
* Path to SAS macro file;					%LET macroloc=  	 C:\ECPE;
* Path to import/export files from; 		%LET filesave=  	 C:\ECPE; 

* Name prefix for files to be created;		%LET filename = 	ECPE_saturated;
* Name of SAS library files are stored in;	%LET saslibname=	work;
* Name of SAS dataset for Q matrix;			%LET Qname= 		ECPEq;
* Name of SAS dataset with original data;	%LET dataname=		ECPEdata_saturated;

* Name of person ID variable (required);	%LET IDname=		ID;
* Item stem in Q matrix (cant be "item");	%LET itemstem=		x;
* List of items to be modeled;				%LET itemlist=		x1-x28;
* Total number of items;					%LET numitem= 		28;
* Variable for order of item model;			%LET ordervar=  	itemorder; 
* Max order of interaction in item model;   %LET maxitemorder= 	2;

* Attribute stem in Q matrix;				%LET attstem=		attribute;
* Number of categories for attributes;      %LET attcat =       2; * currently only set to 2;
* Total number of attributes;				%LET numatt= 		3;
* Number of total classes (2^A);			%LET numclass=		8;
* Use structural model(0=N,1=Y); 			%LET structon=  	0;
* Order of interaction in structural model; %LET structorder= 	2;

* Loosen congergence criteria (0=N,1=Y)?;	%LET loosen=		0;
* Number of processors available for Mplus;	%LET processors=	4;
*******************************************************************************;

* Import Q-matrix into SAS;
DATA &saslibname..&Qname.;
INPUT &itemstem. attribute1-attribute3 &ordervar.;
DATALINES;
1	1	1	0	2
2	0	1	0	1
3	1	0	1	2
4	0	0	1	1
5	0	0	1	1
6	0	0	1	1
7	1	0	1	2
8	0	1	0	1
9	0	0	1	1
10	1	0	0	1
11	1	0	1	2
12	1	0	1	2
13	1	0	0	1
14	1	0	0	1
15	0	0	1	1
16	1	0	1	2
17	0	1	1	2
18	0	0	1	1
19	0	0	1	1
20	1	0	1	2
21	1	0	1	2
22	0	0	1	1
23	0	1	0	1
24	0	1	0	1
25	1	0	0	1
26	0	0	1	1
27	1	0	0	1
28	0	0	1	1
; RUN;

* Import original data into SAS dataset;
DATA &saslibname..&dataname.; 
    INFILE "&filesave.\ecpedata.dat" TRUNCOVER;
	INPUT &IDname. &itemlist.;
RUN;

	
******************************************************************************
**********       			MACRO EXECUTION:		   	   			  ********
**********			NO CHANGES TO MAKE IN THIS SECTION		          ********
******************************************************************************;

* Call master file;
	%INCLUDE "&macroloc.\LCDM_Mplus2.sas";

* Generate Mplus .dat file and input script;
	%CreateMplusInput(
	filesave=&filesave., filename=&filename., saslibname=&saslibname., 
	Qname=&Qname., dataname=&dataname., IDname=&IDname., itemstem=&itemstem.,  
	itemlist=&itemlist., numitem=&numitem., ordervar=&ordervar., maxitemorder=&maxitemorder., 
	attstem=&attstem., attcat=&attcat., numatt=&numatt., numclass=&numclass.,  
	structon=&structon., structorder=&structorder., loosen=&loosen., processors=&processors.);

* Running Mplus  - program to call, input file, output file;
	X CALL "C:\Program Files\Mplus\mplus.exe" 
			"&filesave.\&filename..inp" 
			"&filesave.\&filename..out"; 

* Calling in master program to import Mplus output;
	%ImportMplusOutput;

*** The following datasets have now been created in the work library:
ReadMplus = all Mplus output text in one variable (one row per line)
Respondents = original item responses, prob of class membership, most likely class membership
ClassCounts = model-estimated class counts and proportions
Thresholds = Estimates, SE, Z-score, and p-values per item per threshold
Itemparms = LCDM Intercept (order=0), main effects (order=1), and interactions (order=2+) per item;

*** If using a structural model, these datasets are also created:
ClassMeans = Mean estimates, SE, Z-score, and p-values per class
StrucParms = Intercept (order=0), main effects (order=1), and interactions (order=2+) per attribute;

* Run this to output these tables of results as SAS datasets and as an excel workbook;	
	%WriteData;
