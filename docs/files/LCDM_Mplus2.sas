%MACRO CreateMplusInput(filesave, filename, saslibname, Qname, dataname, 
IDname, itemstem, itemlist, numitem, ordervar, maxitemorder, attstem, attcat, numatt, 
numclass, structon, structorder, loosen, processors);

/* Activate to send log to external file if needed;
PROC PRINTTO log = "&filesave.\&filename._log.txt" NEW; run; */

* Import original data into work library, save as .dat file for Mplus;
DATA _NULL_;
	SET &saslibname..&dataname.;
	FILE "&filesave.\&dataname..dat";
	PUT &IDname. &itemlist.;
RUN;

******************************************************************************
**********       		  Data Manipulation: 						  ********
**********		Getting Q Matrix and Attribute Patterns	 			  ********
******************************************************************************;

* Sorting by original item number;
PROC SORT DATA=&saslibname..&Qname.; BY &itemstem.; RUN;
DATA Qmatrix; SET &saslibname..&Qname.;
	* Renaming attributes to common format;
	ARRAY old(&numatt.) &attstem.1-&attstem.&numatt.;
	ARRAY new(&numatt.) itematt1-itematt&numatt.;
		DO i=1 TO &numatt.; new(i)=old(i); END;
	* Renaming items to common sequential item order;
	item= _N_;
	DROP i &attstem.1-&attstem.&numatt.;
RUN;

* Saving model order per item as macro variables;
%MACRO SaveItemOrder;
DATA _NULL_;
	* Datasets to draw values from;
	%LET dataQmatrix=%SYSFUNC(OPEN(work.Qmatrix,i));
		%DO i=1 %TO &numitem.;
			* Assign new macro variables as global;
			%GLOBAL itemorder&i.;
			* Gets entire row of data;
			%LET row=%SYSFUNC(FETCHOBS(&dataQmatrix.,&i.));
			* Grabs specific value for threshold from that row;
			%LET item&i.=%SYSFUNC(GETVARN(&dataQmatrix.,
					     %SYSFUNC(VARNUM(&dataQmatrix.,&ordervar.))));
			* Transferring value to global macro, and checking via log;
			%LET itemorder&i. = &&item&i.; 
		    %PUT itemorder&i. &&itemorder&i.; 
		%END;
	* Close dataset used in program;
	%LET dataQmatrix=%SYSFUNC(CLOSE(&dataQmatrix.));
RUN;
%MEND SaveItemOrder;
%SaveItemOrder;

/**** Steps in creating a class to pattern table:
* MACRO Classpattern:
* Creates an initial pattern of as many 0s and 1s in a row as needed
* Copies that pattern into separate datasets as needed
* Concatentates (then deletes) the separate datasets
* Result is one dataset per attribute with enough rows for all classes
* MACRO MergeIt:
* Merges (then deletes) datasets from eaach attribute together; */
%MACRO ClassPattern;
	%LET divisor = 2;
	%LET repeat = 1;
	%LET totclass = &numclass.;
	%DO a=1 %TO &numatt.;
		%DO r=1 %TO %EVAL(&repeat.);
			DATA att&a.&r.;
				DO c=1 TO &totclass.;
					IF c LE (&totclass./&divisor.) THEN classatt&a.=0; ELSE classatt&a.=1;
					OUTPUT;
				END; DROP c;
			run;
			%LET repeat = %EVAL(&repeat.+1);
		%END;
		%IF &a.=1 %THEN %LET order = 1; %ELSE %LET order=%EVAL(&order.*2);
		DATA att&a.; SET %DO loop=1 %TO &order.; Att&a.&loop. %END; ; run;
		PROC DATASETS LIB=WORK NOLIST; DELETE %DO loop=1 %TO &order.; Att&a.&loop. %END; ; 
		RUN; QUIT;
		%LET totclass = %EVAL(&totclass./2);
	%END;	
%MEND ClassPattern;
%ClassPattern;
%MACRO MergeIt;
DATA classpattern; RETAIN class; MERGE
	%DO a=1 %TO &numatt.; Att&a. %END; ;
	class = _N_; RUN;
PROC DATASETS LIB=WORK NOLIST; DELETE %DO a=1 %TO &numatt.; Att&a. %END; ; 
RUN; QUIT;
%MEND MergeIt;
%MergeIt;

* Creating item*attribute kernal table from classpattern matrix and Q matrix;
* Result is one row per class per item;
DATA kernel; 
	DO class=1 TO &numclass.;
		DO item=1 TO &numitem.;
			class=class; item=item; OUTPUT;
		END;
	END;
RUN;
DATA kernel; MERGE kernel classpattern; BY class; RUN;
PROC SORT DATA=kernel; BY item class; RUN;
DATA kernel; MERGE kernel Qmatrix; BY item; RUN;

* Scoring by class and item to create kernel threshold values;
* Result is set of variables called "scoreatt" to be activated
  if the item requires the attribute AND the class has it;
DATA kernel; SET kernel;
	ARRAY aclass(&numatt.) classatt1-classatt&numatt.;
	ARRAY aitem(&numatt.)  itematt1-itematt&numatt.;
	ARRAY ascore(&numatt.) scoreatt1-scoreatt&numatt.;
	DO i=1 TO &numatt.;
		IF aclass(i)=1 AND aitem(i)=1 THEN ascore(i)=1; ELSE ascore(i)=0;
	END; DROP i;
RUN;
* Creating single variables that hold the "scoreatt" and attribute patterns;
DATA kernel; LENGTH scorepattern $10 attpattern $&numatt.; SET kernel; 
	scorepattern = CAT(OF item scoreatt1-scoreatt&numatt.);
	attpattern = CAT(OF scoreatt1-scoreatt&numatt.); 
RUN;
PROC SORT DATA=kernel; BY item scorepattern; RUN;
* Creates count variable for whether each "scoreatt" is unique
  to an item to create a new threshold index as needed;
DATA kernel; SET kernel; BY item scorepattern; RETAIN thresh;
		IF FIRST.Item THEN thresh=1;
		IF scorepattern=LAG1(scorepattern) THEN thresh=thresh; 
		ELSE thresh=thresh+1;
RUN;
* The threshold count is one unit too high, so this FIXES it;
DATA kernel; RETAIN class item scorepattern thresh; SET kernel; 
	thresh=thresh-1;
RUN;
PROC SORT DATA=kernel; BY item class; RUN;

******************************************************************************
**********       		  Data Manipulation: 						  ********
**********			Writing Equations out of Patterns		 		  ********
******************************************************************************;

* Creating variables needed for parameters of the LCDM equations;
%MACRO Equations;
DATA kernel;
	* Listing all character variables eventually created;
	LENGTH  NewStruc NewStrucAdd NewParm NewParmAdd 
			Main OrderMain Order12way Order22way
			Order13way Order23way Order33way
			OrderMainMax Order12wayMax Order22wayMax
			Order13wayMax Order23wayMax Order33wayMax $500; 
	SET kernel;
	*** STRUCTURAL MODEL PARAMETERS;
	* Dummy variable to use as placeholder;
	BeginS = .;
	* Main effects;
		%DO a=1 %TO &numatt.;
			IF classatt&a.=1 THEN G_1&a.= "G_1&a.   "; 
				ELSE G_1&a. = G_1&a.;
		%END;
	* Two-way interactions;
	%IF &structorder.>1 %THEN %DO;
			%DO Num1=1 %TO &numatt.; 
				%DO Num2=2 %TO &numatt.;
					%LET check=%EVAL(&Num1.-&Num2.);
					%IF &check.<0 %THEN %DO;
					IF classatt&Num1.=1 AND classatt&Num2.=1 
						THEN G_2&Num1.&Num2.= "G_2&Num1.&Num2.    ";
						ELSE G_2&Num1.&Num2.=  G_2&Num1.&Num2.;
					%END;
				%END; 
			%END; 
	%END;
	* Three-way interactions;  
	%IF &structorder.>2 %THEN %DO;
		%DO Num1=1 %TO &numatt.; 
			%DO Num2=2 %TO &numatt.;
				%DO Num3=3 %TO &numatt.;
					%LET check1=%EVAL(&Num1.-&Num2.);
					%LET check2=%EVAL(&Num2.-&Num3.);
					%IF &check1.<0 AND &check2.<0 %THEN %DO;
					IF classatt&Num1.=1 AND classatt&Num2.=1 AND classatt&Num3.=1
						THEN G_3&Num1.&Num2.&Num3.= "G_3&Num1.&Num2.&Num3.      ";
						ELSE G_3&Num1.&Num2.&Num3.=  G_3&Num1.&Num2.&Num3.;
					%END;
				%END;
			%END; 
		%END; 
	%END;
	* Four-way interactions; 
	%IF &structorder.>3 %THEN %DO;
		%DO Num1=1 %TO &numatt.; 
			%DO Num2=2 %TO &numatt.;
				%DO Num3=3 %TO &numatt.;
					%DO Num4=4 %TO &numatt.;
						%LET check1=%EVAL(&Num1.-&Num2.);
						%LET check2=%EVAL(&Num2.-&Num3.);
						%LET check3=%EVAL(&Num3.-&Num4.);
						%IF &check1.<0 AND &check2.<0 AND &check3.<0 %THEN %DO;
						IF classatt&Num1.=1 AND classatt&Num2.=1 AND classatt&Num3.=1
									AND classatt&Num4.=1
							THEN G_4&Num1.&Num2.&Num3.&Num4.= "G_4&Num1.&Num2.&Num3.&Num4.          ";
							ELSE G_4&Num1.&Num2.&Num3.&Num4.=  G_4&Num1.&Num2.&Num3.&Num4.;
						%END;
					%END;
				%END;
			%END; 
		%END; 
	%END;
	* Five-way interactions;
 	%IF &structorder.>4 %THEN %DO;
		%DO Num1=1 %TO &numatt.; 
			%DO Num2=2 %TO &numatt.;
				%DO Num3=3 %TO &numatt.;
					%DO Num4=4 %TO &numatt.;
						%DO Num5=5 %TO &numatt.;
							%LET check1=%EVAL(&Num1.-&Num2.);
							%LET check2=%EVAL(&Num2.-&Num3.);
							%LET check3=%EVAL(&Num3.-&Num4.);
							%LET check4=%EVAL(&Num4.-&Num5.);
							%IF &check1.<0 AND &check2.<0 AND &check3.<0 AND &check4.<0 %THEN %DO;
							IF classatt&Num1.=1 AND classatt&Num2.=1 AND classatt&Num3.=1
										AND classatt&Num4.=1 AND classatt&Num5.=1
								THEN G_5&Num1.&Num2.&Num3.&Num4.&Num5.= "G_5&Num1.&Num2.&Num3.&Num4.&Num5.          ";
								ELSE G_5&Num1.&Num2.&Num3.&Num4.&Num5.=  G_5&Num1.&Num2.&Num3.&Num4.&Num5.;
							%END;
						%END;
					%END;
				%END;
			%END; 
		%END; 
	%END;
	* Last variable to use as dummy place;
		EndS = .;
	* Creating string variables to use in NEW equations;
		CALL CATX(" ", NewStruc, OF BeginS--EndS);
			NewStruc = TRANWRD(NewStruc, "." , " " );
			NewStruc = STRIP(NewStruc);
		NewStrucAdd = TRANWRD(STRIP(NewStruc), " " , "+" );
	* Creating index for longest possible number of structural parms;
		LongSParm = LENGTH(STRIP(NewStruc));
run;
* Saving longest parm as new variable for class reference statement;
PROC SORT DATA=kernel; BY item DESCENDING LongSParm; RUN;
DATA SaveSLong; SET kernel; BY item;
	IF FIRST.item THEN NewStrucMax=STRIP(NewStruc); ELSE DELETE; RUN;
* Merge longest parm back into data;
DATA kernel; LENGTH NewStrucMax NewStrucMaxAdd $200; MERGE SaveSLong kernel; BY item; 
	NewStrucMaxAdd = TRANWRD(STRIP(NewStrucMax), " " , "+" ); RUN;
* Clearing extra datasets;
PROC DATASETS LIB=WORK NOLIST; DELETE SaveSLong; RUN; QUIT;

PROC SORT DATA=kernel; BY item thresh; RUN;
* Re-opening data;
DATA kernel; SET kernel;
	*** ITEM PARAMETERS;
	* Dummy variable to use as placeholder;
	BeginI=.;
	* Intercepts per item (1);
		%DO i=1 %TO &numitem.; 
			IF item=&i. THEN LInt = "L&i._0    "; 
				ELSE LInt = LInt;
		%END;
	* Main effects per item (up to # attributes);
		%DO Num1=1 %TO &numatt.;
			%DO i=1 %TO &numitem.;
				IF item=&i. AND scoreatt&Num1.=1 THEN L&i._1&Num1. = "L&i._1&Num1.    "; 
					ELSE L&i._1&Num1. = L&i._1&Num1.;
			%END;
		%END;
	* Two-way interactions per item; 
		%DO i=1 %TO &numitem.;
			%DO Num1=1 %TO &numatt.; 
				%DO Num2=2 %TO &numatt.;
					%LET check=%EVAL(&Num1.-&Num2.);
					%IF &check.<0 AND &&itemorder&i.>1 %THEN %DO;
					IF item=&i. AND scoreatt&Num1.=1 AND scoreatt&Num2.=1 
						THEN L&i._2&Num1.&Num2.= "L&i._2&Num1.&Num2.     ";
						ELSE L&i._2&Num1.&Num2.=  L&i._2&Num1.&Num2.;
					%END;
				%END; 
			%END; 
		%END;
	* Three-way interactions per item;
		%DO i=1 %TO &numitem.; 
			%DO Num1=1 %TO &numatt.; 
				%DO Num2=2 %TO &numatt.;
					%DO Num3=3 %TO &numatt.;
						%LET check1=%EVAL(&Num1.-&Num2.);
						%LET check2=%EVAL(&Num2.-&Num3.);
						%IF &check1.<0 AND &check2.<0 AND &&itemorder&i.>2 %THEN %DO;
						IF item=&i. AND scoreatt&Num1.=1 AND scoreatt&Num2.=1 AND scoreatt&Num3.=1
							THEN L&i._3&Num1.&Num2.&Num3.= "L&i._3&Num1.&Num2.&Num3.       ";
							ELSE L&i._3&Num1.&Num2.&Num3.=  L&i._3&Num1.&Num2.&Num3.;
						%END;
					%END;
				%END; 
			%END; 
		%END;
	* Four-way interactions per item; 
		%DO i=1 %TO &numitem.; 
			%DO Num1=1 %TO &numatt.; 
				%DO Num2=2 %TO &numatt.;
					%DO Num3=3 %TO &numatt.;
						%DO Num4=4 %TO &numatt.;
							%LET check1=%EVAL(&Num1.-&Num2.);
							%LET check2=%EVAL(&Num2.-&Num3.);
							%LET check3=%EVAL(&Num3.-&Num4.);
							%IF &check1.<0 AND &check2.<0 AND &check3.<0 AND &&itemorder&i.>3 %THEN %DO;
							IF item=&i. AND scoreatt&Num1.=1 AND scoreatt&Num2.=1 AND scoreatt&Num3.=1
										AND scoreatt&Num4.=1
								THEN L&i._4&Num1.&Num2.&Num3.&Num4.= "L&i._4&Num1.&Num2.&Num3.&Num4.        ";
								ELSE L&i._4&Num1.&Num2.&Num3.&Num4.=  L&i._4&Num1.&Num2.&Num3.&Num4.;
							%END;
						%END;
					%END;
				%END; 
			%END; 
		%END;
	* Five-way interactions per item; 
		%DO i=1 %TO &numitem.;
			%DO Num1=1 %TO &numatt.; 
				%DO Num2=2 %TO &numatt.;
					%DO Num3=3 %TO &numatt.;
						%DO Num4=4 %TO &numatt.;
							%DO Num5=5 %TO &numatt.;
								%LET check1=%EVAL(&Num1.-&Num2.);
								%LET check2=%EVAL(&Num2.-&Num3.);
								%LET check3=%EVAL(&Num3.-&Num4.);
								%LET check4=%EVAL(&Num4.-&Num5.);
								%IF &check1.<0 AND &check2.<0 AND &check3.<0 AND &check4.<0 AND &&itemorder&i.>4 %THEN %DO;
								IF item=&i. AND scoreatt&Num1.=1 AND scoreatt&Num2.=1 AND scoreatt&Num3.=1
											AND scoreatt&Num4.=1 AND scoreatt&Num5.=1
									THEN L&i._5&Num1.&Num2.&Num3.&Num4.&Num5.= "L&i._5&Num1.&Num2.&Num3.&Num4.&Num5.          ";
									ELSE L&i._5&Num1.&Num2.&Num3.&Num4.&Num5.=  L&i._5&Num1.&Num2.&Num3.&Num4.&Num5.;
								%END;
							%END;
						%END;
					%END;
				%END; 
			%END; 
		%END; */
	* Last variable to use as dummy place to end series;
		LastI = .;
	* Creating string variables to use in NEW equations;
		CALL CATX(" ", NewParm, OF BeginI--LastI);
			NewParm = TRANWRD(NewParm, "." , " " );
			NewParm = STRIP(NewParm);
		NewParmAdd = TRANWRD(STRIP(NewParm), " " , "+" );
	* Creating index for longest possible number of parms;
		LongIParm = LENGTH(STRIP(NewParm));
RUN;
* Saving longest parm as new variable for NEW statement;
PROC SORT DATA=kernel; BY item DESCENDING LongIParm; RUN;
DATA SaveILong; SET kernel; BY item; 
	IF FIRST.item THEN NewParmMax=NewParm; ELSE DELETE; RUN;
* Merge longest parm back into data;
DATA kernel; LENGTH NewParmMax $100; MERGE SaveILong kernel; BY item; RUN;
* Clearing extra datasets;
PROC DATASETS LIB=WORK NOLIST; DELETE SaveILong; RUN; QUIT;

* Creating string variables for ordering constraints;
%LET numatt1 = %EVAL(&numatt.-1);
%LET numatt2 = %EVAL(&numatt.-2);
DATA kernel; SET kernel;
	* Main effects;
		CALL CATX(" ", Main, OF L1_11--L&numitem._1&numatt.);
		* Combining into one variable;
		OrderMain = CATT(TRANWRD(STRIP(Main), " " , ">0; " ), ">0; " );
		* Creating index for longest possible number of parms;
		LongOrderMain = LENGTH(STRIP(OrderMain));
	* 1 Two-way interactions - need two sets of loops for both possible constraints;
		* Creating dummy variable to begin first two-way interaction series;
		BeginO12way=.;
		%DO i=1 %TO &numitem.;
			%DO Num1=1 %TO &numatt.; 
				%DO Num2=2 %TO &numatt.;
					%LET check1=%EVAL(&Num1.-&Num2.);
					%IF &check1.<0 AND &&itemorder&i.>1 %THEN %DO;
						IF item=&i. AND scoreatt&Num1.=1 AND scoreatt&Num2.=1 THEN DO;
							O1L&i._2&Num1.&Num2. = "L&i._2&Num1.&Num2.>-L&i._1&Num1.;         " ;
						END;
					%END;
				%END;
			%END;
		%END;
		* Creating dummy variable to end first two-way interaction series;
		EndO12way=.;
	* 2 Two-way interactions - need two sets of loops for both possible constraints;
		* Creating dummy variable to begin second two-way interaction series;
		BeginO22way=.;
		%DO i=1 %TO &numitem.;
			%DO Num1=1 %TO &numatt.; 
				%DO Num2=2 %TO &numatt.;
					%LET check1=%EVAL(&Num1.-&Num2.);
					%IF &check1.<0 AND &&itemorder&i.>1 %THEN %DO;
						IF item=&i. AND scoreatt&Num1.=1 AND scoreatt&Num2.=1 THEN DO;
							O2L&i._2&Num1.&Num2. = "L&i._2&Num1.&Num2.>-L&i._1&Num2.;         " ;
						END;
					%END;
				%END;
			%END;
		%END;
		* Creating dummy variable to end second two-way interaction series;
		EndO22way=.;
		%IF &maxitemorder.>1 %THEN %DO;
			* Combining into one variable;
			Order12way = CATS(OF BeginO12way--EndO12way);
			Order12way = TRANWRD(Order12way, ".", " ");
			Order22way = CATS(OF BeginO22way--EndO22way);
			Order22way = TRANWRD(Order22way, ".", " ");
			* Creating index for longest possible number of parms;
			LongOrder12way = LENGTH(STRIP(Order12way));
			LongOrder22way = LENGTH(STRIP(Order22way));
		%END;
	* 1 Three-way interactions - need three sets of loops for all three sets of constraints;
		* Creating dummy variable to begin first three-way interaction series;
		BeginO13way=.;
		%DO i=1 %TO &numitem.;
			%DO Num1=1 %TO &numatt.; 
				%DO Num2=2 %TO &numatt.;
					%DO Num3=3 %TO &numatt.;
						%LET check1=%EVAL(&Num1.-&Num2.);
						%LET check2=%EVAL(&Num2.-&Num3.);
						%IF &check1.<0 AND &check2.<0 AND &&itemorder&i.>2 %THEN %DO;
							IF item=&i. AND scoreatt&Num1.=1 AND scoreatt&Num2.=1 AND scoreatt&Num3.=1 THEN DO;
								O1L&i._3&Num1.&Num2.&Num3. = "L&i._3&Num1.&Num2.&Num3.>-(L&i._2&Num2.&Num3.+L&i._2&Num1.&Num3.+L&i._1&Num3.);	     	  " ;
							END;
						%END;
					%END;
				%END;
			%END;
		%END;
		* Creating dummy variable to end first three-way interaction series;
		EndO13way=.;
	* 2 Three-way interactions - need three sets of loops for all three sets of constraints;
		* Creating dummy variable to begin second three-way interaction series;
		BeginO23way=.;
		%DO i=1 %TO &numitem.;
			%DO Num1=1 %TO &numatt.; 
				%DO Num2=2 %TO &numatt.;
					%DO Num3=3 %TO &numatt.;
						%LET check1=%EVAL(&Num1.-&Num2.);
						%LET check2=%EVAL(&Num2.-&Num3.);
						%IF &check1.<0 AND &check2.<0 AND &&itemorder&i.>2 %THEN %DO;
							IF item=&i. AND scoreatt&Num1.=1 AND scoreatt&Num2.=1 AND scoreatt&Num3.=1 THEN DO;
								O2L&i._3&Num1.&Num2.&Num3. = "L&i._3&Num1.&Num2.&Num3.>-(L&i._2&Num2.&Num3.+L&i._2&Num1.&Num2.+L&i._1&Num2.);	     	  " ;
							END;
						%END;
					%END;
				%END;
			%END;
		%END;
		* Creating dummy variable to end second three-way interaction series;
		EndO23way=.;
	* 3 Three-way interactions - need three sets of loops for all three sets of constraints;
		* Creating dummy variable to begin third three-way interaction series;
		BeginO33way=.;
		%DO i=1 %TO &numitem.;
			%DO Num1=1 %TO &numatt.; 
				%DO Num2=2 %TO &numatt.;
					%DO Num3=3 %TO &numatt.;
						%LET check1=%EVAL(&Num1.-&Num2.);
						%LET check2=%EVAL(&Num2.-&Num3.);
						%IF &check1.<0 AND &check2.<0 AND &&itemorder&i.>2 %THEN %DO;
							IF item=&i. AND scoreatt&Num1.=1 AND scoreatt&Num2.=1 AND scoreatt&Num3.=1 THEN DO;
								O3L&i._3&Num1.&Num2.&Num3. = "L&i._3&Num1.&Num2.&Num3.>-(L&i._2&Num1.&Num3.+L&i._2&Num1.&Num2.+L&i._1&Num1.);	     	  " ;
							END;
						%END;
					%END;
				%END;
			%END;
		%END;
		* Creating dummy variable to end third three-way interaction series;
		EndO33way=.;
		%IF &maxitemorder.>2 %THEN %DO;
			* Combining into one variable;
			Order13way = CATS(OF BeginO13way--EndO13way);
			Order13way = TRANWRD(Order13way, '.' , ' ' );
			Order23way = CATS(OF BeginO23way--EndO23way);
			Order23way = TRANWRD(Order23way, '.' , ' ' );
			Order33way = CATS(OF BeginO33way--EndO33way);
			Order33way = TRANWRD(Order33way, '.' , ' ' );
			* Creating index for longest possible number of parms;
			LongOrder13way = LENGTH(STRIP(Order13way));
			LongOrder23way = LENGTH(STRIP(Order23way));
			LongOrder33way = LENGTH(STRIP(Order33way));
		%END; 
RUN;
*** Saving longest parm as new variable for ordering statements;
* Main effects;
	PROC SORT DATA=kernel; BY item DESCENDING LongOrderMain; RUN;
	DATA SaveMain; SET kernel; BY item; 
		* Save max value;
		IF FIRST.item THEN OrderMainMax=STRIP(OrderMain); ELSE DELETE; 
	RUN;
	* Merge longest parm back into data;
	DATA kernel; MERGE SaveMain kernel(DROP=OrderMainMax); BY item; RUN;
	* Clearing extra datasets;
	PROC DATASETS LIB=WORK NOLIST; DELETE SaveMain; RUN; QUIT;
* Two-way interactions;
%IF &maxitemorder. >1 %THEN %DO; 
* Do for each set; 
* Set 1;
	PROC SORT DATA=kernel; BY item DESCENDING LongOrder12way; RUN;
	DATA Save2way; SET kernel; BY item;
		* Save max value; 
		IF FIRST.item THEN Order12wayMax=STRIP(Order12way); ELSE DELETE; 
		* Add blanks between semi-colons;
		Order12wayMax = TRANWRD(Order12wayMax, ';' , '; ' ); 
	RUN;
	* Merge longest parm back into data;
	DATA kernel; MERGE Save2way kernel(DROP=Order12wayMax); BY item; RUN;
* Set 2;
	PROC SORT DATA=kernel; BY item DESCENDING LongOrder22way; RUN;
	DATA Save2way; SET kernel; BY item;
		* Save max value; 
		IF FIRST.item THEN Order22wayMax=STRIP(Order22way); ELSE DELETE; 
		* Add blanks between semi-colons;
		Order22wayMax = TRANWRD(Order22wayMax, ';' , '; ' ); 
	RUN;
	* Merge longest parm back into data;
	DATA kernel; MERGE Save2way kernel(DROP=Order22wayMax); BY item; RUN;
	* Clearing extra datasets;
	PROC DATASETS LIB=WORK NOLIST; DELETE Save2way; RUN; QUIT;
%END;
* Three-way interactions;
%IF &maxitemorder. >2 %THEN %DO; 
* Do for each set; 
* Set 1;
	PROC SORT DATA=kernel; BY item DESCENDING LongOrder13way; RUN;
	DATA Save3way; SET kernel; BY item;
		* Save max value; 
		IF FIRST.item THEN Order13wayMax=STRIP(Order13way); ELSE DELETE; 
		* Add blanks between semi-colons;
		Order13wayMax = TRANWRD(Order13wayMax, ';' , '; ' ); RUN;
	* Merge longest parm back into data;
	DATA kernel; MERGE Save3way kernel(DROP=Order13wayMax); BY item; RUN;
* Set 2;
	PROC SORT DATA=kernel; BY item DESCENDING LongOrder23way; RUN;
	DATA Save3way; SET kernel; BY item;
		* Save max value; 
		IF FIRST.item THEN Order23wayMax=STRIP(Order23way); ELSE DELETE; 
		* Add blanks between semi-colons;
		Order23wayMax = TRANWRD(Order23wayMax, ';' , '; ' ); RUN;
	* Merge longest parm back into data;
	DATA kernel; MERGE Save3way kernel(DROP=Order23wayMax); BY item; RUN;
* Set 2;
	PROC SORT DATA=kernel; BY item DESCENDING LongOrder33way; RUN;
	DATA Save3way; SET kernel; BY item;
		* Save max value; 
		IF FIRST.item THEN Order33wayMax=STRIP(Order33way); ELSE DELETE; 
		* Add blanks between semi-colons;
		Order33wayMax = TRANWRD(Order33wayMax, ';' , '; ' ); RUN;
	* Merge longest parm back into data;
	DATA kernel; MERGE Save3way kernel(DROP=Order33wayMax); BY item; RUN;
* Clearing extra datasets;
	PROC DATASETS LIB=WORK NOLIST; DELETE Save3way; RUN; QUIT;
%END;
PROC SORT DATA=kernel; BY item class; RUN;
%MEND Equations;

%Equations;

*** Creating datasets for use in writing script;
* Full item by class data for writing class-based commands;
DATA itembyclass; SET kernel; RUN;
PROC SORT DATA=itembyclass; BY class item; RUN;
* Class-level data for writing structural model;
DATA classlevel; SET kernel; WHERE item=1; RUN;
PROC SORT DATA=classlevel; BY class; RUN;
* Item by threshold data for writing item-based commands;
DATA itembythresh; SET kernel; RUN;
PROC SORT NODUPKEY DATA=itembythresh; BY item thresh; RUN;
DATA itembythresh; SET itembythresh; DROP class classatt1-classatt&numatt.; RUN;
	

******************************************************************************
**********       		  Writing Mplus Script 						  ********
******************************************************************************;

* Writing Mplus script for DCM;
%MACRO WriteMplus;
DATA _NULL_;
* Datasets to draw values from;
%LET dataitembyclass=%SYSFUNC(OPEN(work.itembyclass,i));
%LET dataitembythresh=%SYSFUNC(OPEN(work.itembythresh,i));
%LET dataclasslevel=%SYSFUNC(OPEN(work.classlevel,i));

* Name of Mplus file to be written with returns and breaks;
    FILE "&filesave.\&filename..inp" PRINT;
* Mplus TITLE command;
	PUT @1 "TITLE:  ! Section that appears in header of output file";
	PUT	@5 "DCM for &dataname. with &numatt. attributes and &structorder.-order structural model,";
	PUT	@5 "&numitem. items, and maximum &maxitemorder.-order item model,";
	PUT	@5 "Saturated structural model (Mplus default)." /;
* Mplus DATA command;
	PUT	@1 "DATA:  ! Location of free format data file" 
	/	@5 "FILE = &filesave.\&dataname..dat;"
	/ 	;
* Mplus VARIABLE command;
	PUT @1 "VARIABLE:"
	/	@5 "NAMES = &IDname. mitem1-mitem&numitem.;" @40 "! List of variables in data file" 
	/	@5 "USEVARIABLE = mitem1-mitem&numitem.;"    @40 "! Variables to be analyzed" 
	/	@5 "CATEGORICAL = mitem1-mitem&numitem.;"    @40 "! Binary outcomes" 
	/	@5 "CLASSES = c(&numclass.);" @40 "! Number of possible attribute patterns (2^A)" 
	/	@5 "IDvariable = &IDname.;"   @40 "! Person ID variable to save respondent data"
	/	;
* Mplus ANALYSIS command;
	PUT @1 "ANALYSIS:" 
	/	@5 "TYPE = MIXTURE;" @40 "! Estimates latent classes" 
	/	@5 "STARTS = 0;"     @40 "! Turn off multiple random start feature" 
	/	@5 "PROCESSORS = &processors.;" @40 "! Number of processors available"
%IF &loosen.=1 %THEN %DO;
	/	@5 "mconvergence=  .1;"  @40 "! Less stringent convergence criteria"
    /	@5 "muconvergence= .1;"
    /	@5 "convergence=   .1;"
    /	@5 "logcriterion=  .1;"
	/	@5 "rlogcriterion= .1;"
    /	@5 "h1convergence= .1;"
    /	@5 "mcconvergence= .1;"
    /	@5 "miterations=   10000;"
	/	@5 "LOGHIGH=	   10;"
	/	@5 "LOGLOW=       -10;"
%END;
	/	;
* Mplus MODEL command;
	PUT @1 "MODEL:" 
	/	;

* Code for structural model if needed;
	%IF &structon.=1 %THEN %DO;
		PUT	@1 "%" "OVERALL" "%";
	* List means per class; 
		%DO c=1 %TO %EVAL(&numclass.-1);
			PUT @1 "[c#&c.] (m&c.);  ! Latent variable mean for class &c.";
		%END; PUT / ;
	%END;

* List class-specific models; 
	%LET counter=1;
	%DO c=1 %TO &numclass.;
		PUT @1 "%" "c#&c." "%"  "  ! Model for Class &c.";
			%DO i=1 %TO &numitem.;
				* Gets entire row of data;
				%LET row=%SYSFUNC(FETCHOBS(&dataitembyclass.,&counter.));
				* Grabs specific value for threshold from that row;
				%LET tnum=%SYSFUNC(GETVARN(&dataitembyclass,
						  %SYSFUNC(VARNUM(&dataitembyclass,thresh))));
				PUT @5 "[mitem&i.$1] (T&i._&tnum.);" @30 "! Item &i. Thresh &tnum.";
				%LET counter=%EVAL(&counter.+1); 
			%END;
		PUT / ;
	%END; 

* Mplus MODEL CONSTRAINT and NEW commands; 
	PUT @1 "MODEL CONSTRAINT:  ! Used to define LCDM parameters"
	/	@1 "! Mplus uses P(X=0) rather than P(X=1) so multiply by -1"
	/	;

* Code for structural model if needed;
	%IF &structon.=1 %THEN %DO;
	PUT	@1 "! STRUCTURAL MODEL";
	* List class-specific means; 
		%LET c=1;
			* Gets entire row of data;
			%LET row=%SYSFUNC(FETCHOBS(&dataclasslevel.,&c.));
			* Grabs specific value for threshold from that row;
			%LET newstrucline=%SYSFUNC(GETVARC(&dataclasslevel.,
					  %SYSFUNC(VARNUM(&dataclasslevel.,NewStrucMax))));
			PUT @1 "NEW(G_0 &newstrucline.);";
		%DO c=1 %TO %EVAL(&numclass.-1);
			* Gets entire row of data;
			%LET row2=%SYSFUNC(FETCHOBS(&dataclasslevel.,&c.));
			* Grabs specific value for threshold from that row;
			%LET strucadd=%SYSFUNC(GETVARC(&dataclasslevel.,
					  %SYSFUNC(VARNUM(&dataclasslevel,NewStrucAdd))));
			%LET strucmaxadd=%SYSFUNC(GETVARC(&dataclasslevel.,
					  %SYSFUNC(VARNUM(&dataclasslevel,NewStrucMaxAdd))));
			PUT @1 "m&c.= &strucadd. - (&strucmaxadd.);";
		%END; 
			PUT @1 "G_0= - (&strucmaxadd.);" /;
	%END;

* NEW commands per item; 
	%LET counter=1;
	%LET counterd= %EVAL(&counter.-1);
	* Have to do it manually the first time through;
	%LET i=1;
		* Gets entire row of current data;
		%LET nowrow=%SYSFUNC(FETCHOBS(&dataitembythresh.,&counter.));
		* Grabs specific value for variable from that row;
		%LET newline=%SYSFUNC(GETVARC(&dataitembythresh,
					 %SYSFUNC(VARNUM(&dataitembythresh,NewParmMax))));
		PUT @1 "! Item &i: Define LCDM parameters present for item &i.";
		PUT @1 "NEW(&newline.);";
		*** Equations per item;
		* Have to do first threshold manually;
		%LET t=1;
			* Gets entire row of current data;
			%LET nowrow=%SYSFUNC(FETCHOBS(&dataitembythresh.,&counter.));
			* Grabs specific value for variable from that row;
			%LET threshline=%SYSFUNC(GETVARC(&dataitembythresh,
						 %SYSFUNC(VARNUM(&dataitembythresh,NewParmAdd))));
			PUT @1  "T&i._&t.=-(&threshline.);" @72 "! Item &i. Thresh &t.";
			%LET counter=%EVAL(&counter.+1); %LET counterd= %EVAL(&counter.-1);
		* Iteratively through rest of thresholds;
		%DO t=2 %TO &numclass.;
			* Gets entire row of current data;
			%LET nowrow=%SYSFUNC(FETCHOBS(&dataitembythresh.,&counter.));
			* Grabs specific value for variable from that row;
			%LET nowitem=%SYSFUNC(GETVARN(&dataitembythresh,
						 %SYSFUNC(VARNUM(&dataitembythresh,item))));
			* Grabs specific value for variable from that row;
			%LET threshline=%SYSFUNC(GETVARC(&dataitembythresh,
						    %SYSFUNC(VARNUM(&dataitembythresh,NewParmAdd))));
			* Gets entire row of previous data;
			%LET pastrow=%SYSFUNC(FETCHOBS(&dataitembythresh.,&counterd.));
			* Grabs specific value for variable from that row;
			%LET pastitem=%SYSFUNC(GETVARN(&dataitembythresh,
						  %SYSFUNC(VARNUM(&dataitembythresh,item))));
				%IF &nowitem.=&pastitem. %THEN %DO;
					PUT @1  "T&i._&t.=-(&threshline.);" @72 "! Item &i. Thresh &t.";
					%LET counter=%EVAL(&counter.+1); %LET counterd= %EVAL(&counter.-1);
				%END;
		%END;
		*** Ordering constraints per item;
		* Gets entire row of current data;
		%LET nowrow=%SYSFUNC(FETCHOBS(&dataitembythresh.,&counterd.));
		* Grabs specific value for variable from that row;
		%LET listmain=%SYSFUNC(GETVARC(&dataitembythresh,
					  %SYSFUNC(VARNUM(&dataitembythresh,OrderMainMax))));
		* Grabs specific value for variable from that row;
		%IF &&itemorder&i.>1 %THEN %DO;
		%LET list12way=%SYSFUNC(GETVARC(&dataitembythresh,
					   %SYSFUNC(VARNUM(&dataitembythresh,Order12wayMax))));
		%LET list22way=%SYSFUNC(GETVARC(&dataitembythresh,
					   %SYSFUNC(VARNUM(&dataitembythresh,Order22wayMax))));
		%END;
		%IF &&itemorder&i.>2 %THEN %DO;
		%LET list13way=%SYSFUNC(GETVARC(&dataitembythresh,
					   %SYSFUNC(VARNUM(&dataitembythresh,Order13wayMax))));
		%LET list23way=%SYSFUNC(GETVARC(&dataitembythresh,
					   %SYSFUNC(VARNUM(&dataitembythresh,Order23wayMax))));
		%LET list33way=%SYSFUNC(GETVARC(&dataitembythresh,
					   %SYSFUNC(VARNUM(&dataitembythresh,Order33wayMax))));
		%END;
		* Need to check to see if all terms exist, otherwise get "" in input;
		PUT @1 "! Main effect order constraints"
		/   @1 "&listmain.";
		%IF &&itemorder&i.>1 AND %LENGTH(&list12way.)>0 %THEN %DO;
		PUT @1 "! Two-way interaction order constraints"
		/   @1 "&list12way."
		/   @1 "&list22way."; %END;
		%IF &&itemorder&i.>2 AND %LENGTH(&list13way.)>0 %THEN %DO;
		PUT @1 "! Three-way interaction order constraints"
		/   @1 "&list13way."
		/   @1 "&list23way."
		/   @1 "&list33way."; %END;
		PUT / ;

	* Iterially through rest of items, checking for new thresholds to write;
	%DO i=2 %TO &numitem.;
		* Gets entire row of current data;
		%LET nowrow=%SYSFUNC(FETCHOBS(&dataitembythresh.,&counter.));
		* Grabs specific value for variable from that row;
		%LET newline=%SYSFUNC(GETVARC(&dataitembythresh,
					 %SYSFUNC(VARNUM(&dataitembythresh,NewParmMax))));
		PUT @1 "! Item &i: Define LCDM parameters present for item &i.";
		PUT @1 "NEW(&newline.);";
		*** Equations per item;
		* Have to do first threshold manually;
		%LET t=1;
			* Gets entire row of current data;
			%LET nowrow=%SYSFUNC(FETCHOBS(&dataitembythresh.,&counter.));
			* Grabs specific value for variable from that row;
			%LET threshline=%SYSFUNC(GETVARC(&dataitembythresh,
						 %SYSFUNC(VARNUM(&dataitembythresh,NewParmAdd))));
			PUT @1 "T&i._&t.=-(&threshline.);" @72 "! Item &i. Thresh &t.";
			%LET counter=%EVAL(&counter.+1); %LET counterd= %EVAL(&counter.-1);
		* Iteratively through rest of items;
		%DO t=2 %TO &numclass.;
			* Gets entire row of current data;
			%LET nowrow=%SYSFUNC(FETCHOBS(&dataitembythresh.,&counter.));
			* Grabs specific value for variable from that row;
			%LET nowitem=%SYSFUNC(GETVARN(&dataitembythresh,
						 %SYSFUNC(VARNUM(&dataitembythresh,item))));
			* Grabs specific value for variable from that row;
			%LET threshline=%SYSFUNC(GETVARC(&dataitembythresh,
						 %SYSFUNC(VARNUM(&dataitembythresh,NewParmAdd))));
			* Gets entire row of previous data;
			%LET pastrow=%SYSFUNC(FETCHOBS(&dataitembythresh.,&counterd.));
			* Grabs specific value for variable from that row;
			%LET pastitem=%SYSFUNC(GETVARN(&dataitembythresh,
						  %SYSFUNC(VARNUM(&dataitembythresh,item))));
				%IF &nowitem.=&pastitem. %THEN %DO;
					PUT @1 "T&i._&t.=-(&threshline.);" @72 "! Item &i. Thresh &t.";
					%LET counter=%EVAL(&counter.+1); %LET counterd= %EVAL(&counter.-1);
				%END;
		%END;
		*** Ordering constraints per item;
		* Gets entire row of current data;
		%LET nowrow=%SYSFUNC(FETCHOBS(&dataitembythresh.,&counterd.));
		* Grabs specific value for variable from that row;
		%LET listmain=%SYSFUNC(GETVARC(&dataitembythresh,
					  %SYSFUNC(VARNUM(&dataitembythresh,OrderMainMax))));
		* Grabs specific value for variable from that row;
		%IF &&itemorder&i.>1 %THEN %DO;
		%LET list12way=%SYSFUNC(GETVARC(&dataitembythresh,
					   %SYSFUNC(VARNUM(&dataitembythresh,Order12wayMax))));
		%LET list22way=%SYSFUNC(GETVARC(&dataitembythresh,
					   %SYSFUNC(VARNUM(&dataitembythresh,Order22wayMax))));
		%END;
		%IF &&itemorder&i.>2 %THEN %DO;
		%LET list13way=%SYSFUNC(GETVARC(&dataitembythresh,
					   %SYSFUNC(VARNUM(&dataitembythresh,Order13wayMax))));
		%LET list23way=%SYSFUNC(GETVARC(&dataitembythresh,
					   %SYSFUNC(VARNUM(&dataitembythresh,Order23wayMax))));
		%LET list33way=%SYSFUNC(GETVARC(&dataitembythresh,
					   %SYSFUNC(VARNUM(&dataitembythresh,Order33wayMax))));
		%END;
		* Need to check to see if all terms exist, otherwise get "" in input;
		PUT @1 "! Main effect order constraints"
		/   @1 "&listmain.";
		%IF &&itemorder&i.>1 AND %LENGTH(&list12way.)>0 %THEN %DO;
		PUT @1 "! Two-way interaction order constraints"
		/   @1 "&list12way."
		/   @1 "&list22way."; %END;
		%IF &&itemorder&i.>2 AND %LENGTH(&list13way.)>0 %THEN %DO;
		PUT @1 "! Three-way interaction order constraints"
		/   @1 "&list13way."
		/   @1 "&list23way."
		/   @1 "&list33way."; %END;
		PUT / ;
	%END; 

* Mplus OUTPUT command;
	PUT @1 "OUTPUT:"
	/	@5 "TECH10;  ! Request additional model fit statistics"
	/	;

* Mplus SAVEDATA command;
	PUT @1 "SAVEDATA: ! Format, name of posterior probabilities of class membership file"
	/	@5 "FORMAT = F10.5;"
	/	@5 "FILE = &filesave.\respondents_&dataname..dat;"
	/	@5 "SAVE = CPROBABILITIES;"
	/	;
run;

* Close datasets used in program;
%LET dataitembyclass=%SYSFUNC(CLOSE(&dataitembyclass.));
%LET dataitembythresh=%SYSFUNC(CLOSE(&dataitembythresh.));
%LET dataclasslevel=%SYSFUNC(CLOSE(&dataclasslevel.));

%MEND WriteMplus;

%WriteMplus; 

/* Close log files being printed to if had activated;
	PROC PRINTTO; run; */

%MEND CreateMplusInput;


%MACRO ImportMplusOutput;

* Reading in person-level output;
DATA respondents; RETAIN &IDname.;
	INFILE "&filesave.\respondents_&dataname..dat" DLM="     " ;
	INPUT &itemlist. &IDname. cprob1-cprob&numclass. class;
RUN;

DATA ClassPattern2; SET ClassPattern; DROP class; RUN;

* This is to compute the probability of mastery for each examinee;
PROC IML;
	USE ClassPattern2;
	READ ALL VAR _ALL_ INTO pattmat;

	USE Respondents;
	READ ALL VAR {&IDname.} INTO idvec;
	nexam=NROW(idvec);

	cproblist=CONCAT("cprob",CHAR(1,1,0));
	DO i=2 to &numclass.;
		         IF i < 10  THEN cproblist=cproblist||CONCAT("cprob",CHAR(i,1,0));
			ELSE IF I < 100 THEN cproblist=cproblist||CONCAT("cprob",CHAR(i,2,0));
				            ELSE cproblist=cproblist||CONCAT("cprob",CHAR(i,3,0));
	END;

	USE respondents;
	READ ALL VAR cproblist INTO probmat;

	attprob=probmat*pattmat;

	attprob=idvec||attprob;

	cname={"&IDname."};
	
	DO i=1 TO &numatt.;
		IF i < 10 THEN;	cname=cname//CONCAT("prob_&attstem.",CHAR(i,1,0));
			   ELSE 	cname=cname//CONCAT("prob_&attstem.",CHAR(i,2,0));
	END;

	CREATE attprob FROM attprob [colname=cname];
	APPEND FROM attprob;

QUIT;

* Trimming up respondents file -- cutting out items and adding attribute probabilities;
PROC SORT DATA=respondents; BY &IDname.; RUN;
PROC SORT DATA=attprob; BY &IDname.; RUN;
DATA respondents (DROP=&itemlist.); MERGE respondents attprob; BY &IDname.; RUN; */


* Reading ALL output into a one-variable file;
DATA ReadMplus;
	INFILE "&filesave.\&filename..out" DSD TRUNCOVER LRECL=9000;
	INPUT line $200.;
run;

* Do if used a structural model;
%IF &structon.=1 %THEN %DO;
	* Getting class means;
	DATA ClassMeans; SET ReadMplus;
		* Selecting only class means;
			WHERE INDEX(line, "#");
		* Deleting lines from model input;
			IF INDEX(line, "!") THEN DELETE;
		* Extracting values into separate variables;
			class = 	INPUT(SUBSTR(line, ANYDIGIT(line)), 8.0);
			CMeanEst = INPUT(SCAN(line,2, " "), 10.3);
			CMeanSE =  INPUT(SCAN(line,3, " "), 10.3);
			CMeanZ =   INPUT(SCAN(line,4, " "), 10.3);
			CMeanP =   INPUT(SCAN(line,5, " "), 10.3);
		DROP line;
	RUN;
	* Getting NEW parameters for structural model;
	DATA StrucParms; SET ReadMplus;
		* Selecting only class means;
			WHERE INDEX(line, "G") AND INDEX(line, "_");
		* Deleting lines from model input;
			IF INDEX(line, ")") THEN DELETE;
		* Extracting values into separate variables;
			parm = SCAN(line,1, " ");
			Order =    INPUT(SUBSTR(line, ANYPUNCT(line)+1, 1), 8.0);
			Atts = 	   INPUT(SUBSTR(line, ANYPUNCT(line)+2), 8.0);
			IF Order = 0 THEN Atts = 0;
			StrucEst = INPUT(SCAN(line,2, " "), 10.3);
			StrucSE =  INPUT(SCAN(line,3, " "), 10.3);
			StrucZ =   INPUT(SCAN(line,4, " "), 10.3);
			StrucP =   INPUT(SCAN(line,5, " "), 10.3);
		DROP line;
	RUN;
%END;

* Getting estimated model class counts and proportions;
DATA ClassCounts; SET ReadMplus;
	temp = SUBSTR(line, 1, FIND(line," "));
	%DO i=1 %TO &numclass.;
		IF temp=&i. THEN DO;
			class = INPUT(temp, 3.0);
			estcount = INPUT(SCAN(line,2, " "), 10.3);
			estprop =  INPUT(SCAN(line,3, " "), 10.3);
		END;
	%END;
	IF NMISS(estcount,estprop)>0 THEN DELETE;
	KEEP class estcount estprop;
RUN; 
DATA ClassCounts; SET ClassCounts;
	IF _N_ > &numclass. THEN DELETE;
RUN;

* Merge class attribute patterns into class counts;
DATA ClassCounts; MERGE ClassCounts ClassPattern; BY class; RUN;

* Calculate reliabilities;
PROC IML;
	USE ClassCounts;
	READ ALL VAR {estprop} INTO cprob;
	
	cname="prob_&attstem.1";
	DO i=2 TO &numatt.;
		IF i < 10 THEN;	cname=cname||CONCAT("prob_&attstem.",CHAR(i,1,0));
			   ELSE 	cname=cname||CONCAT("prob_&attstem.",CHAR(i,2,0));
	END;
   
	USE Respondents;
	READ ALL VAR cname INTO rprob;
	nexam=NROW(rprob);

	USE ClassPattern2;
	READ ALL VAR _ALL_ INTO pattmat;
	
	attprob=cprob`*pattmat;  

	reliability=J(1,&numatt.,0);
	DO i=1 TO nexam;
		DO j=1 to &numatt.;
			vartrue=attprob[1,j]*(1-attprob[1,j]);
			varerror=rprob[i,j]*(1-rprob[i,j]);
			reliability[1,j]=reliability[1,j]+(vartrue/(vartrue+varerror));
		END;
	END;
	reliability=reliability/nexam;
	
	attnumber=J(&numatt.,1,0);
	DO i=1 TO &numatt.;
		attnumber[i,1]=i;
	END;
	outmat=attnumber||attprob`||reliability`;
	
	CREATE AttProbRel FROM outmat [colname = {Attribute Probability Reliability}];
	APPEND FROM outmat;

QUIT;

* Needed for parsing item number out of threshold label;
      %IF &numitem.<10  %THEN %LET digit=1; 
%ELSE %IF &numitem.<100 %THEN %LET digit=2;
%ELSE                         %LET digit=3;

* Getting thresholds per class per item;
DATA Thresholds; SET ReadMplus;
	* Selecting only thresholds;
		WHERE INDEX(line, "$");
	* Deleting lines with comments (from model input);
		IF INDEX(line, "!") THEN DELETE;
	* Extracting values into separate variables;
		Item = 		INPUT(SUBSTR(line, 6, &digit.), 8.0);
		Thresh = 	INPUT(SUBSTR(line, ANYPUNCT(line)+1, FIND(line," ")), 8.0);
		ThreshEst = INPUT(SCAN(line,2, " "), 10.3);
		ThreshSE =  INPUT(SCAN(line,3, " "), 10.3);
		ThreshZ =   INPUT(SCAN(line,4, " "), 10.3);
		ThreshP =   INPUT(SCAN(line,5, " "), 10.3);
RUN;
DATA Thresholds; RETAIN class item thresh; SET Thresholds;
	* Adding index for class;
	IF _N_=1 THEN class = 1; 
	IF item > LAG1(item) THEN class=class; ELSE class=class+1;
	DROP line;
RUN;
* Merge class pattern and item attribute pattern info;
DATA Thresholds; MERGE Thresholds ItembyClass 
	(KEEP= class item &itemstem. itematt1-itematt&numatt. classatt1-classatt&numatt.);
	BY class item; RUN;
DATA Thresholds; RETAIN class item &itemstem. thresh; SET Thresholds; RUN;

* Getting NEW parameters for item model;
DATA ItemParms; SET ReadMplus;
	* Selecting only class means;
	WHERE INDEX(line, "L") AND INDEX(line, "_");
	* Deleting lines from model input;
	IF INDEX(line, ")") OR INDEX(line, ">") OR INDEX(line, "FILE") THEN DELETE;
	parm = SCAN(line,1, " ");
	Item = 		INPUT(SUBSTR(line, ANYDIGIT(line), ANYPUNCT(line)-2), 8.0);
	Order = 	INPUT(SUBSTR(line, ANYPUNCT(line)+1, 1), 8.0);
	Atts = 		INPUT(SUBSTR(line, ANYPUNCT(line)+2), 8.0);
	IF Order = 0 THEN Atts = 0;
	ItemEst = INPUT(SCAN(line,2, " "), 10.3);
	ItemSE =  INPUT(SCAN(line,3, " "), 10.3);
	ItemZ =   INPUT(SCAN(line,4, " "), 10.3);
	ItemP =   INPUT(SCAN(line,5, " "), 10.3);
	DROP line;
RUN;
* Merge back original item labels;
DATA ItemParms; MERGE ItemParms Qmatrix; BY Item; RUN;
DATA ItemParms; RETAIN parm item &itemstem. order; SET ItemParms; RUN;


* Clearing extra datasets;
PROC DATASETS LIB=WORK NOLIST; 
	DELETE Classlevel ClassPattern ClassPattern2 ItembyClass ItembyThresh kernel 
		   AttProb; 
RUN; QUIT;

%MEND ImportMplusOutput;



%MACRO WriteData;
* Run to export summary datasets to permanent library named above;
	DATA folder.&filename._Qmatrix;     SET Qmatrix;     RUN;
	DATA folder.&filename._ReadMplus;   SET ReadMplus;   RUN;
	DATA folder.&filename._Respondents; SET Respondents; RUN;
	DATA folder.&filename._ClassCounts; SET ClassCounts; RUN;
	DATA folder.&filename._Thresholds;  SET Thresholds;  RUN;
	DATA folder.&filename._ItemParms;   SET ItemParms;   RUN;
	DATA folder.&filename._AttProbRel;  SET AttProbRel;  RUN;

* Run to export to an excel workbook (one dataset per worksheet);
	PROC EXPORT DATA=Qmatrix     OUTFILE= "&filesave.\&filename..xls" DBMS=EXCEL REPLACE; SHEET="Qmatrix";     RUN;
	PROC EXPORT DATA=ReadMplus   OUTFILE= "&filesave.\&filename..xls" DBMS=EXCEL REPLACE; SHEET="ReadMplus";   RUN;
	PROC EXPORT DATA=Respondents OUTFILE= "&filesave.\&filename..xls" DBMS=EXCEL REPLACE; SHEET="Respondents"; RUN;
	PROC EXPORT DATA=ClassCounts OUTFILE= "&filesave.\&filename..xls" DBMS=EXCEL REPLACE; SHEET="ClassCounts"; RUN;
	PROC EXPORT DATA=Thresholds  OUTFILE= "&filesave.\&filename..xls" DBMS=EXCEL REPLACE; SHEET="Thresholds";  RUN;
	PROC EXPORT DATA=ItemParms   OUTFILE= "&filesave.\&filename..xls" DBMS=EXCEL REPLACE; SHEET="ItemParms";   RUN;
	PROC EXPORT DATA=AttProbRel  OUTFILE= "&filesave.\&filename..xls" DBMS=EXCEL REPLACE; SHEET="Attribute Reliability";   RUN;

* Additional output if using a structural model;
%IF &structon.=1 %THEN %DO;
* Run to export summary datasets to permanent library named above;
	DATA folder.&filename._Classmeans;  SET ClassMeans;  RUN;
	DATA folder.&filename._StrucParms;  SET StrucParms;  RUN;
* Run to export to an excel workbook (one dataset per worksheet);
	PROC EXPORT DATA=ClassMeans  OUTFILE= "&filesave.\&filename..xls" DBMS=EXCEL REPLACE; SHEET="ClassMeans";  RUN;
	PROC EXPORT DATA=StrucParms  OUTFILE= "&filesave.\&filename..xls" DBMS=EXCEL REPLACE; SHEET="StrucParms";  RUN;
%END;

%MEND WriteData; 
