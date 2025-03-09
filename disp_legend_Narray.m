% iLt:index of legend term
iLt=0;
Lts=Narray;%可调项
% Lts=Vs;%可调项
t1='';t2='';t3='';t4='';t5='';t6='';t7='';t8='';t9='';t10='';
for iLt=1:length(Lts)
    t0='';
     t0=strcat(t0,'L=');
%     t0=strcat(t0,'V=');
    t0=strcat(t0,num2str(Lts(iLt)));
    if(iLt==1)t1=t0;end
    if(iLt==2)t2=t0;end
    if(iLt==3)t3=t0;end
    if(iLt==4)t4=t0;end
    if(iLt==5)t5=t0;end
    if(iLt==6)t6=t0;end
    if(iLt==7)t7=t0;end
    if(iLt==8)t8=t0;end
    if(iLt==9)t9=t0;end
    if(iLt==10)t10=t0;end
end
legendtext=legend(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10);
set(legendtext,'box','off');
set(legendtext,'FontSize',15);
