% PHplot.m %
clear;
n=400;
k=1.1;
l=1;
load phases;
load inputdistr;
a=phases(:,1)';
p=length(a);
T=phases(:, 2:p+1);

t = -T*ones(size(a'));
vv = -a*inv(T)*ones(size(a'));
var = 2*a*inv(T*T)*ones(size(a'))-vv*vv;
std=sqrt(var);
truncpoint = vv+4*sqrt(var);
dt = truncpoint/n;
for i=1:n+1, 
  y(i)=dt*(i-1);
  F(i)=1-a*expm(T*y(i))*ones(size(a'));
  S(i)=1-F(i);
  f(i)=a*expm(T*y(i))*t;
  r(i)=f(i)/S(i);
end
for i=1:length(inputdistr(:,1)),
  x1(i)=inputdistr(i,1);
  F1(i)=inputdistr(i,2);
  S1(i)=1-F1(i);
  f1(i)=inputdistr(i,3);
  r1(i)=inputdistr(i,4);
end
while l
  disp(' ');
  disp('1. Mean and standard-deviation');
  disp('2. Display survival function');
  disp('3. Display distribution function');
  disp('4. Display density');
  disp('5. Display failure rate');
  disp('6. Print graph');
  disp('7. Save graph on file "PHgraph"'); 
  disp('8. Load new estimates of pi and T'); 
  disp('9. Quit');
    
  q=input('Select (1-9): ');
  
  if q==1
    disp(' ');
    disp('The fitted PH-distribution has mean ='), disp(vv); 
    disp('and standard-deviation ='), disp(std); 
  end
  if q==2
    plot(y,S,'--',x1,S1,'-');
    title('Survival function');
    xlabel('-  input, - - fitted PH');
    ylabel(' ');
    set(gca,'xlim',[0,y(n)],'ylim',[0,1]);
  end
  if q==3
    plot(y,F,'--',x1,F1,'-');
    title('Distribution function');
    xlabel('-  input, - - fitted PH');
    ylabel(' ');
    set(gca,'xlim',[0,y(n)],'ylim',[0,1]);
  end
  if q==4
    plot(y,f,'--',x1,f1,'-');
    title('Density');
    xlabel('- - fitted PH');
    ylabel(' ');
    set(gca,'xlim',[0,y(n)],'ylim',[0,max(max(f), max(f1))*k]);
 end
  if q==5
    plot(y,r,'--',x1,r1,'-');
    title('Failure rate');
    xlabel('- - fitted PH');
    ylabel(' ');
    set(gca,'xlim',[0,y(n)],'ylim',[0,max(max(r),max(r1))*k]);
  end
  if q==6
    print;
  end
  if q==7
    print PHgraph -dps;
  end
  if q==8
    clear a T t;
    clear y F S f r;
    load phases;
    a=phases(:,1)';
    p=length(a);
    T=phases(:, 2:p+1);
    
    t = -T*ones(size(a'));
    vv = -a*inv(T)*ones(size(a'));
    var = 2*a*inv(T*T)*ones(size(a'))-vv*vv;
    std=sqrt(var);
    truncpoint = vv+4*sqrt(var);
    dt = truncpoint/n;
    for i=1:n+1, 
      y(i)=dt*(i-1);
      F(i)=1-a*expm(T*y(i))*ones(size(a'));
      S(i)=1-F(i);
      f(i)=a*expm(T*y(i))*t;
      r(i)=f(i)/S(i);
    end  
  end
  if q==9
    l=0;
  end
end
