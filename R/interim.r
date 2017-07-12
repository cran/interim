#' Scheduling interim analyses in clinical trials
#'
#' It is often discussed during the planning of a clinical trial whether an interim analysis
#' is beneficial. The time point of the interim  analysis and the end of the clinical trial
#' are crucial for the decision. Both  depend on the recruitment of patients and on the length
#' of the treatment phase. The package \code{interim} allows the instantaneous simulation and
#' plotting of both the recruitment and treatment phase. Based on these simulations, the
#' timing of interim analyses can be assessed or different recruitment scenarios
#' can be compared.
#'
#' @details
#' There are three main functions in this package:
#' \itemize{
#' \item \code{\link{recruitment}}
#' \item \code{\link{treatment}}
#' \item \code{\link{trialCourse}}
#' }
#'
#' The function \code{recruitment} generally is the starting point. It simulates screening
#' and enrollment based on screening characteristics, like e.g. number of available centers
#' and patients, screen failure rate, etc.
#'
#' The function \code{treatment} simulates the treatment period based on a given recruitment
#' scenario as simulated by \code{recruitment}.
#'
#' The function \code{trialCourse} plots displays of enrollment and treatment simulations.
#' Two plots are provided; the first one displays center openings and the second one
#' displays patient screening, enrollment and treatment.
#'
#' In addition to the main functions, the package comprises a number of auxilliary functions
#' helping to derive or convert parameters required for the three main functions,
#' as well as to derive and plot information on timing of interim analyses. The following
#' auxilliary functions are available:
#' \itemize{
#' \item \code{\link{capacity}}
#' \item \code{\link{convertedRate}}
#' \item \code{\link{trialWeek}}
#' \item \code{\link{cross}}
#' }


"_PACKAGE"



  Umap <- function(f,x,...) unlist(Map(f,x,...))

  yield <- function(t,oc,cc,wc) list(weeksOfTrial=t,openCenters=oc,closedCenters=cc,centerWeeks=wc)

  ### open centers
   f <- function(x,x0,b) if (x <= x0) b*x else b*x0

   ### closed centers
   g <- function(x,x0,b) if (x <= x0) 0 else b*(x-x0)


   ### lower left triangle
   L <- function(x,b) b/2*x^2

   l <- function(b,L) sqrt(2*L/b)


   ### parallelogram in the middle
   M <- function(x,xl,h) h*(x-xl)

   m <- function(xl,h,M) M/h+xl


   ### upper right triangle
   R <- function(x,xl,xr,b) (x-xl)*b*(xr-xl)-b/2*(x-xl)^2

   r <- function(xl,xr,b,R) xr-sqrt((xr-xl)^2-2*R/b)


   w1 <- function(cw,z){
      te=l(cw,z);
      t=unique(sort(c(0:floor(te),te)))
      oc=Umap(f,t,Inf,cw)
      cc=Umap(g,t,Inf,cw)
      wc=Umap(L,t,cw)
      yield(t,oc,cc,wc)
   }


   w2 <- function(nc,cw,z){
      tz=nc/cw
      h=L(tz,cw)
      if (z <= h)
         w1(cw,z)
      else {
         te=m(tz,nc,z-h)
         t=unique(sort(c(0:floor(te),te,tz)))
         oc=Umap(f,t,tz,cw)
         cc=Umap(g,t,Inf,cw)
         wL=Umap(L,t[t <= tz],cw)
         wM=Umap(M,t[t > tz],tz,cw*tz)
         wc=c(wL,h+wM)
         yield(t,oc,cc,wc)
      }
    }


   w3 <- function(cw,ns,sw,z){
      ta=ns/sw
      h=L(ta,cw)
      if (z <= h)
         w1(cw,z)
      else {
         te=m(ta,cw*ta,z-h)
         t=unique(sort(c(0:floor(te),te,ta)))
         oc=Umap(f,t,Inf,cw)
         cc=Umap(g,t,ta,cw)
         wL=Umap(L,t[t <= ta],cw)
         wM=Umap(M,t[t > ta],ta,cw*ta)
         wc=c(wL,h+wM)
         yield(t,oc,cc,wc)
      }
   }



   w4 <- function(nc,cw,ns,sw,z){
      tz=nc/cw; ta=ns/sw
      t1=min(tz,ta); t2=max(tz,ta);     t3=tz+ta
      h1=L(t1,cw);   h2=M(t2,t1,cw*t1); h3=R(t3,t2,t3,cw)

      if (z <= h1) ### te in lower left triangle
         w1(cw,z)
      else
      if (z <= h1+h2) ### te in parallelogram
         if (tz <=  ta)
            w2(cw,nc,z)
         else
            w3(cw,ns,sw,z)
      else
      if (z <= h1+h2+h3) { ### te in upper right triangle
         te=r(t2,t3,cw,z-h1-h2)
         t=unique(sort(c(0:floor(te),te,t2,t1)))
         oc=Umap(f,t,tz,cw)
         cc=Umap(g,t,ta,cw)
         wL=Umap(L,t[t <= t1],cw)
         wR=Umap(R,t[t > t2],t2,t3,cw)
         if (t1 < t2) { ### size of parallelogram greater than zero
            if (tz < ta)
               hP=nc    ### steep parallelogram
            else
               hP=cw*ta ### slant parallelogram
            wM=Umap(M,t[(t1 < t) & (t <= t2)],t1,hP)
            wc=c(wL,h1+wM,h1+h2+wR)
         } else { ### size of parallelogram equal to zero
            wc=c(wL,h1+wR) ###  this clause is needed because of vector length
         }
         yield(t,oc,cc,wc)
      } else
      yield(NA,NA,NA,NA)
   }


   bccBlue=rgb(0,144,197,maxColorValue=255)
   bccGreen=rgb(107,194,0,maxColorValue=255)
   bccGrey=rgb(103,103,103,maxColorValue=255)

   coCol="#5ca754" # green
   ccCol="#cb4b41" # red
   scCol="#165b97" # darkblue
   enCol="#597dae" # midblue
   t1Col="#81a5c9" # lightblue
   t2Col="#81a5c9" # lightblue



##############
### EXPORT ###
##############



#' Scheduling interim analyses in clinical trials
#'
#' Function \code{recruitment} simulates the recruitment of patients in clinical trials.
#' @param nc maximum number of centers to be opened or \code{Inf}.
#' @param ns maximum number of patients to be screened within each center or \code{Inf}.
#' @param cw number of center openings per week.
#' @param sw number of screened patients per week within each center.
#' @param sf screening failure rate.
#' @param tb time between screening and enrollment/randomization in weeks.
#' @param en number of patients to be enrolled.
#'
#' @details
#' Function \code{recruitment} simulates the recruitment progress for a required
#' number of enrolled patients in clinical trials based on the expected number
#' of centers to be opened per week and the expected number of patients being recruited per site and week.
#' The function assumes that centers are being opened at a constant rate per week (\code{cw})
#' and patient per center are screened at a constant rate per week (\code{sw}).
#'
#' The function can handle recruitment limits by limiting the total number of centers (\code{nc})
#' and/or the number of patients recruitable per site (\code{ns}).
#'
#' The function discriminates between screening timepoint and enrollment/randomization timepoint by
#' allowing a screening period (\code{tb}) and screen failure rate (\code{sf}) to be specified.
#' If both are zero then patients are directly enrolled at screening.
#'
#' Function \code{recruitment} can handle four different recruitment scenarios.
#'
#' \itemize{
#'    \item Scenario 1. Centers are being opened over the entire trial duration,
#' i.e. no limit of centers (\code{nc=Inf}) and are kept open during the complete trial duration,
#' i.e. no limit of patients per center (\code{ns=Inf}).
#'    \item Scenario 2. Only a limited number of centers can be opened (\code{nc=} a positive number) and are
#' kept open during the complete trial duration (\code{ns=Inf}).
#'    \item Scenario 3. Centers are being opened over the entire trial duration (\code{nc=Ind}) and are
#' only have a limited capacity for patient recruitment (\code{ns=} a positive number).
#'    \item Scenario 4. Only a limited number of centers can be opened (\code{nc=} a positive number) and are
#' only have a limited capacity for patient recruitment (\code{ns=} a positive number).
#' }
#'
#' Under scenario 4 only a limited number of patients can be recruited. The auxilliary function \code{capacity} can
#' be used to derive the maximum number of patients that can be enrolled under this scenario.
#'
#' Results of \code{recruitment} are required as input for function \code{\link{treatment}} to derive
#' the time when treatment of patients is finished.
#'
#' @return
#' \code{recruitment} returns a list of vectors with the following components:
#' \itemize{
#' \item \code{weeksOfTrial} a vector counting the trial week for recruitment(from 0 (start of site openings) to the required number of weeks for recruitment)
#' \item \code{openCenters} a vector with the (cumulative) number of opened centers per trial week
#' \item \code{closedCenters} a vector with the (cumulative) number of closed centers per trial week
#' \item \code{centerWeeks} a vector with the (cumulative) number of opened center weeks per trial week
#' \item \code{screenings} a vector with the (cumulative) number of screened patients per trial week
#' \item \code{weeksOfEnrollment} a vector counting the weeks of enrollment (with start of site openings as reference start)
#' \item \code{enrollments} a vector with the (cumulative) number of enrolled/randomized patients per week of enrollment
#' }
#' @export recruitment
#'
#' @seealso
#' \code{\link{treatment}} for simulating the treatment duration for a given recruitment scenario;
#' \code{\link{trialCourse}} for plots of recruitment and treatment scenarios;
#' \code{\link{capacity}} for deriving the maximum number of patients that can be enrolled under scenario 4;
#'
#' @examples
#' x=recruitment(nc=Inf,ns=Inf,cw=4,sw=2,sf=0.3,tb=4,en=400)

   recruitment <- function(nc,ns,cw,sw,sf,tb,en){
      z=en/(1-sf)/sw

      ### model 1
      if (is.infinite(nc) & is.infinite(ns)) h=w1(cw,z) else
      ### model 2
      if (is.finite(nc)   & is.infinite(ns)) h=w2(nc,cw,z) else
      ### model 3
      if (is.infinite(nc) & is.finite(ns))   h=w3(cw,ns,sw,z) else
      ### model 4
      if (is.finite(nc)   & is.finite(ns))   h=w4(nc,cw,ns,sw,z)

      h$screenings=sw*h$centerWeeks
      h$weeksOfEnrollment=h$weeksOfTrial+tb
      h$enrollments=(1-sf)*h$screenings

      h
   }



#' Scheduling interim analyses in clinical trials
#'
#' Function \code{capacity} calculates the maximum number of enrollments for a recruitment in scenario 4.
#' @param nc maximum number of centers to be opened.
#' @param ns maximum number of patients to be screened within each center.
#' @param sf screening failure rate.
#'
#' @details
#' \code{capacity} is an auxilliary function to determine the maximum number of patients that can be enrolled
#' in the scenario where only a limited number of centers is available and each center only has a limited number
#' of patients that can be enrolled.
#'
#' @return
#' The maximum number of patients that can be enrolled.
#'
#' @export capacity
#'
#' @seealso
#' \code{\link{recruitment}} for simulating recruitment scenarios
#'
#' @examples
#' mE=capacity(nc=40,ns=10,sf=0.3)
#' recruitment(nc=40,ns=10,cw=4,sw=2,sf=0.3,tb=4,en=mE)

   capacity <- function(nc,ns,sf){
      (1-sf)*nc*ns
   }



#' Scheduling interim analyses in clinical trials
#'
#' Function \code{convertedRate} converts a drop-out rate from one period to another. If the drop-out rate is defined for
#' w1 weeks \code{convertedRate} yields the drop-out rate for w2 weeks.
#'
#' @param r rate between 0 and 1 (0 < r < 1)
#' @param w1 number of weeks for which r is defined
#' @param w2 number of weeks to which the rate shall be converted
#'
#' @details
#' \code{convertedRate} is an auxilliary function that converts drop-out rates for different time periods.
#'
#' The function can be used in order to convert drop-out rates required for function \code{treatment}. Function
#' \code{treatment} requires the drop-out rate for the respective treatment duration as input. Typically known annual
#' drop-out rates can be converted to the respective rate for the treatment duration accordingly by setting \code{w1}
#' to 52 and \code{w2} to the respective treatment duration.
#'
#' @return
#' The converted drop-out rate.
#'
#' @export convertedRate
#'
#' @seealso
#' \code{\link{treatment}} for simulating the treatment duration for a given recruitment scenario
#'
#' @examples
#' convertedRate(r=0.3,w1=52,w2=26)

   convertedRate <- function(r,w1,w2){
      d=1-(1-r)^(1/w1) ## d is the weekly rate
      1-(1-d)^w2
   }




#' Scheduling interim analyses in clinical trials
#'
#' Function \code{treatment} simulates the treatment phase base on a recruitment scenario simulated by function \code{recruitment}.
#' @param r recruitment scenario calculated with function \code{recruitment}.
#' @param du duration of treatment phase in weeks.
#' @param dr drop-out rate during the treatment phase.
#'
#' @details
#' \code{treatment} simulates the treatment period based on a given recruitment scenario.
#' The function assumes a common fixed treatment length for all subjects (\code{du}).
#' The drop-out rate may be included via \code{dr}. If drop-out rates are available
#' only for different time periods, e.g. annual rates, function \code{\link{convertedRate}} can be used to convert
#' the rate to a drop-out rate for the respective treatment duration.
#'
#' @return
#' \itemize{
#' \item \code{treatment} returns a list of vectors with the following components:
#' \item \code{patients} a vector with the (cumulative) number of patients who finished treatment
#' \item \code{weeksOfTrial} a vector with the corresponding trial week when patients have finished treatment
#' (with start of site openings as reference start)
#' }
#' @export treatment
#'
#' @seealso
#' \code{\link{recruitment}} for simulating recruitment scenarios;
#' \code{\link{trialCourse}} for plots of recruitment and treatment scenarios;
#' \code{\link{convertedRate}} for converting drop-out rates for different time periods.
#'
#' @examples
#' x=recruitment(nc=Inf,ns=Inf,cw=4,sw=2,sf=0.3,tb=4,en=400)
#' y=treatment(r=x,du=26,dr=0.2)
#' z=treatment(r=x,du=52,dr=0.2)
#'
#' x=recruitment(nc=Inf,ns=Inf,cw=4,sw=2,sf=0.3,tb=4,en=400)
#' y=treatment(r=x,du=26,dr=convertedRate(0.3,52,26))
#' z=treatment(r=x,du=52,dr=0.3)

   treatment <- function(r,du,dr){
      h=list()

      h$patients=(1-dr)*r$enrollments
      h$weeksOfTrial=r$weeksOfEnrollment+du

      h
   }




#' Scheduling interim analyses in clinical trials
#'
#' Function \code{trialCourse} plots the results of function \code{recruitment}
#' and function \code{treatment}.
#' @param r recruitment scenario calculated with function \code{recruitment}.
#' @param t1 \emph{optional}. Treatment phase simulation from function \code{treatment}.
#' @param t2 \emph{optional}. Treatment phase simulation from function \code{treatment}.
#' @param lp \emph{optional}. Position of legend, specified by keyword: "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right", or "center".
#'
#' @details
#' Function \code{trialCourse} produces two plots to display results of enrollment
#' and treatment simulations.
#'
#' The first plot displays the cumulative number of centers that have been opened
#' as well as the cumulative number of centers that have been closed, if applicable, per trial week.
#'
#' The second plot displays the number of patients that have been screened and enrolled per trial week.
#' If at least one of the parameters \code{t1} and \code{t2} are not \code{NULL}, then
#' the number of patients finished treatment is also displayed.
#'
#' It is possible to include two different treatment scenarios into one plot. This option may for example
#' be used to assess the timing for specific interim analyses, i.e. \code{t1} is used to assess when the
#' required number of patients for the interim analysis finished treatment while \code{t2} is used to assess
#' when the required number of patients for the final analysis finished treatment.
#'
#' @export trialCourse
#'
#' @import graphics
#' @seealso
#' \code{\link{treatment}} for simulating the treatment duration for a given recruitment scenario;
#' \code{\link{recruitment}} for simulating recruitment scenarios.
#'
#' @examples
#' x=recruitment(nc=Inf,ns=Inf,cw=4,sw=2,sf=0.3,tb=4,en=400)
#' y=treatment(r=x,du=26,dr=convertedRate(0.3,52,26))
#' z=treatment(r=x,du=52,dr=0.3)
#' trialCourse(r=x,t1=y,t2=z)

   trialCourse <- function(r,t1=NULL,t2=NULL,lp="topright"){

      layout(matrix(c(1,2),nrow=1,ncol=2))

      w=union(r$weeksOfTrial,r$weeksOfEnrollment)
      yLabel="Screened and enrolled patients"

      if (!is.null(t1)) {
         w=union(w,t1$weeksOfTrial)
         yLabel="Screened, enrolled, and treated patients"
      }

      if (!is.null(t2)) {
         w=union(w,t2$weeksOfTrial)
         yLabel="Screened, enrolled, and treated patients"
      }

      ### center plot
      y=c(rep(0,length(w)-1),max(r$openCenters))
      plot(w,y,type="n",main="Centers",xlab="Week",ylab="Openings and closings")
      lines(r$weeksOfTrial,r$openCenters,ylim=c(0,max(r$openCenters,r$closedCenter)),type="l",lwd=3,col=coCol)
      lines(r$weeksOfTrial,r$closedCenters,type="l",lwd=3,col=ccCol)

      legend(lp,lwd=3,col=c(coCol,ccCol),legend=c("Center openings","Center closings"))

      ### patient plot
      y=c(rep(0,length(w)-1),max(r$screenings))
      plot(w,y,type="n",main="Patients",xlab="Week",ylab=yLabel)
      lines(r$weeksOfTrial,r$screenings,type="l",lwd=3,col=scCol)
      lines(r$weeksOfEnrollment,r$enrollments,type="l",lwd=3,col=enCol)
      if (!is.null(t1))
         lines(t1$weeksOfTrial,t1$patients,type="l",lty=6,lwd=3,col=t1Col)
      if (!is.null(t2))
         lines(t2$weeksOfTrial,t2$patients,type="l",      lwd=3,col=t2Col)

      if (is.null(t1) & is.null(t2))
         legend(lp,lwd=3,col=c(scCol,enCol),legend=c("Screened patients","Enrolled patients"))
      else
      if (!is.null(t1) & is.null(t2))
      	 legend(lp,lwd=3,lty=c(1,1,6),col=c(scCol,enCol,t1Col),legend=c("Screened patients","Enrolled patients","Treated patients 1"))
      else
      if (is.null(t1) & !is.null(t2))
      	 legend(lp,lwd=3,col=c(scCol,enCol,t2Col),legend=c("Screened patients","Enrolled patients","Treated patients 2"))
      else
         legend(lp,lwd=3,lty=c(1,1,6,1),col=c(scCol,enCol,t1Col,t2Col),legend=c("Screened patients","Enrolled patients", "Treated patients 1","Treated patients 2"))
   }




#' Scheduling interim analyses in clinical trials
#'
#' Function \code{trialWeek} determines the week of the trial in which a certain number \code{p}
#' of patients finished treatment.
#' @param t result of function \code{treatment}.
#' @param p number of patients for which the week shall be determined.
#'
#' @details
#' \code{trialWeek} is an auxilliary function required to assess the timing of interim analyses. It derives
#' the week of trial in which a certain number of patients finished treatment.
#'
#' The output is required for function \code{cross}, which includes the information into an existing Patients diagram.
#'
#' @return
#' The week in which the number of patients is reached.
#'
#' @export trialWeek
#' @seealso
#' \code{\link{cross}} for plotting results of function \code{trialWeek} into an existing Patients diagram.
#'
#' @examples
#' x=recruitment(nc=Inf,ns=Inf,cw=4,sw=2,sf=0.3,tb=4,en=400)
#' y=treatment(r=x,du=26,dr=convertedRate(0.3,52,26))
#' z=treatment(r=x,du=52,dr=0.3)
#' trialCourse(r=x,t1=y,t2=z)
#' trialWeek(t=y,p=100)

   trialWeek <- function(t,p){
      x=t$weeksOfTrial
      y=t$patients

      if (p < y[1])
         NA
      else
         if (p == y[1])
            x[1]
         else {
            n=length(y)
            if (p <= y[n]){
               s=(1:(n-1))[(y[-n] < p) & (p <= y[-1])] # segment
               y1=y[s]; y2=y[s+1]
               x1=x[s]; x2=x[s+1]
               b=(y2-y1)/(x2-x1)
               a=((y1+y2)-b*(x1+x2))/2
               (p-a)/b
            } else
               NA
         }
   }




#' Scheduling interim analyses in clinical trials
#'
#' Function \code{cross} plots two cossing lines into the patients diagram .
#' @param w week where the vertical line is plotted.
#' @param p number of patients where the horizontal line is plotted.
#'
#' @details
#' This function includes a vertical and horizontal line into an existing patient diagram
#' produced by function \code{trialCourse}.
#'
#' The lines are to mark the timepoint \code{w} in weeks at which a required number of
#' patients \code{p} has finished their treatment. The display can be used to assess
#' the scheduling of interim analyses.
#'
#' The auxilliary function \code{\link{trialWeek}} can be used to derive the week of the
#' trial in which the required number of patients has finished the treatment.
#'
#' @export cross
#'
#' @seealso
#' \code{\link{trialCourse}} for plots of recruitment and treatment scenarios;
#' \code{\link{trialWeek}} for deriving the week of a trial at which a certain number of patients finished treatment.
#'
#' @examples
#' x=recruitment(nc=Inf,ns=Inf,cw=4,sw=2,sf=0.3,tb=4,en=400)
#' y=treatment(r=x,du=26,dr=convertedRate(0.3,52,26))
#' z=treatment(r=x,du=52,dr=0.3)
#' trialCourse(r=x,t1=y,t2=z)
#' trialWeek(t=y,p=100)
#' cross(w=trialWeek(t=y,p=100),p=100)

   cross <- function(w,p){
      abline(v=w,col="black")
      abline(h=p,col="black")
   }



###############
### E N D E ###
###############
