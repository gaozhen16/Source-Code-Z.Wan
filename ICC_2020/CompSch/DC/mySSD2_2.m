function [ AbhhHat0,cout] = mySSD2_2( yy,W,N1,N2,L,maxNum )
    cout = 0;
    yytemp = yy;
    ppos = [];
    maxDis = ceil((maxNum-4)/8);
for ss = 1:L
    qq = zeros(N1-1,1);
    for q = 1+maxDis : N1-1-maxDis
        Wq = W(:,(q-1)*N2+1:(q+1)*N2);
        qq(q) = norm(Wq \ yytemp);
    end
    [~,sq] = max(qq);

    pp = zeros(N2-1,1);
    for p = 1+maxDis : N2-1-maxDis
        Wp = W(:,unique([p:N2:p+(N1-1)*N2,p+1:N2:p+1+(N1-1)*N2]));
        pp(p) = norm(Wp \ yytemp);
    end
    [~,sp] = max(pp);
    maxY = sq;
    maxX = sp;
    xCou = 1;
    yCou = 1;
    maxDir = 5;
    maxCou = 0;
    secY = [maxY+yCou maxY      maxY+yCou];
    secX = [maxX      maxX+xCou maxX+xCou];
    secDir = [5 5 6];
    secCou = [yCou xCou xCou];
    AllY = [maxY,secY];
    AllX = [maxX,secX];
    AllDir = [maxDir secDir];
    AllCou = [maxCou secCou];
    forBack = zeros(4,3);
    layer = 2;
    dir = 1;
    cou = -1;

    i = 4;
    for layer = 1:N1*N2
        for ccout = 1:8
            if ccout ==1 AllY = [AllY AllY(ccout  )-yCou*layer];AllX = [AllX AllX(ccout)             ];dir = yCou+3;cou = 0;   end
            if ccout ==2 AllY = [AllY AllY(ccout  )+yCou*layer];AllX = [AllX AllX(ccout)             ];dir = -yCou+3;cou = 0;   end
            if ccout ==3 AllY = [AllY AllY(ccout  )-yCou*layer];AllX = [AllX AllX(ccout)             ];dir = yCou+3;cou = xCou;end
            if ccout ==4 AllY = [AllY AllY(ccout  )+yCou*layer];AllX = [AllX AllX(ccout)             ];dir = -yCou+3;cou = xCou;end
            if ccout ==5 AllY = [AllY AllY(ccout-4)           ];AllX = [AllX AllX(ccout-4)-xCou*layer];dir = xCou+2;cou = 0;   end
            if ccout ==6 AllY = [AllY AllY(ccout-4)           ];AllX = [AllX AllX(ccout-4)-xCou*layer];dir = xCou+2;cou = yCou;end
            if ccout ==7 AllY = [AllY AllY(ccout-4)           ];AllX = [AllX AllX(ccout-4)+xCou*layer];dir = -xCou+2;cou = 0;   end
            if ccout ==8 AllY = [AllY AllY(ccout-4)           ];AllX = [AllX AllX(ccout-4)+xCou*layer];dir = -xCou+2;cou = yCou;end
            AllDir = [AllDir dir];
            AllCou = [AllCou cou];
            i = i+1;
            if i == maxNum break;end
        end
        if i == maxNum break;end
    end

AllPos = AllX + N2 * (AllY - 1);

AbhhHat = W(:,AllPos) \ yytemp;
AbhhHat0 = zeros(N1*N2,1);
AbhhHat0(AllPos) = AbhhHat;

% PlotH(AbhhHat0);

outPos = AllPos(maxNum-8+1:end);
AbhhHatOut = AbhhHat0(outPos);
[~,tempMaxPos0] = max(AbhhHatOut);
[~,tempMinPos0] = min(AbhhHatOut);
tempMaxPosAll = AllPos(find(AbhhHat==AbhhHatOut(tempMaxPos0)));
tempMinPosAll = AllPos(find(AbhhHat==AbhhHatOut(tempMinPos0)));
tempMaxPosPar = find(AbhhHat==AbhhHatOut(tempMaxPos0));
tempMinPosPar = find(AbhhHat==AbhhHatOut(tempMinPos0));
forBack( AllDir(tempMaxPosPar),AllCou(tempMaxPosPar)+2 ) = 1;
forBack( AllDir(tempMinPosPar),AllCou(tempMinPosPar)+2 ) = -1;
while 1
    if AllDir(tempMaxPosPar) ==1 newY = AllY(tempMaxPosPar)  ;newX = AllX(tempMaxPosPar)+1;end
    if AllDir(tempMaxPosPar) ==2 newY = AllY(tempMaxPosPar)+1;newX = AllX(tempMaxPosPar)  ;end
    if AllDir(tempMaxPosPar) ==3 newY = AllY(tempMaxPosPar)  ;newX = AllX(tempMaxPosPar)-1;end
    if AllDir(tempMaxPosPar) ==4 newY = AllY(tempMaxPosPar)-1;newX = AllX(tempMaxPosPar)  ;end
    
    if newX>N2 || newX<1 || newY>N1 || newY<1
        break;
    end
    
    AllY = [AllY newY];
    AllX = [AllX newX];
    AllCou = [AllCou,AllCou(tempMaxPosPar)];
    AllDir = [AllDir AllDir(tempMaxPosPar)];
    AllPos = [AllPos newX + N2 * (newY - 1)];
    outPos(outPos==tempMaxPosAll) = [];
    outPos = [outPos newX + N2 * (newY - 1)];
    
    if AllDir(tempMinPosPar) ==1 newY2 = AllY(tempMinPosPar)  ;newX2 = AllX(tempMinPosPar)-1;end
    if AllDir(tempMinPosPar) ==2 newY2 = AllY(tempMinPosPar)-1;newX2 = AllX(tempMinPosPar)  ;end
    if AllDir(tempMinPosPar) ==3 newY2 = AllY(tempMinPosPar)  ;newX2 = AllX(tempMinPosPar)+1;end
    if AllDir(tempMinPosPar) ==4 newY2 = AllY(tempMinPosPar)+1;newX2 = AllX(tempMinPosPar)  ;end
    if AllDir(tempMinPosPar) ==5 newY2 = []                   ;newX2 = []                   ;end
    if AllDir(tempMinPosPar) ==6 newY2 = []                   ;newX2 = []                   ;end
    AllPos(tempMinPosPar) = [];
    AllX(tempMinPosPar) = [];
    AllY(tempMinPosPar) = [];
    AllDir(tempMinPosPar) = [];
    AllCou(tempMinPosPar) = [];
    outPos(outPos==tempMinPosAll) = [];
    if ~isempty(newX2 + N2 * (newY2 - 1)) 
        if ~(isempty(find(AllPos==newX2 + N2 * (newY2 - 1)))) 
            outPos = [outPos newX2 + N2 * (newY2 - 1)];
        end
    end
    
    cout = cout + 1;
    AbhhHat = W(:,AllPos) \ yytemp;
    AbhhHat0 = zeros(N1*N2,1);
    AbhhHat0(AllPos) = AbhhHat;
    
%     PlotH(AbhhHat0);
    
    AbhhHatOut = AbhhHat0(outPos);
    [~,tempMaxPos0] = max(AbhhHatOut);
    [~,tempMinPos0] = min(AbhhHatOut);
    tempMaxPosAll = AllPos(find(AbhhHat==AbhhHatOut(tempMaxPos0)));
    tempMinPosAll = AllPos(find(AbhhHat==AbhhHatOut(tempMinPos0)));
    tempMaxPosPar = find(AbhhHat==AbhhHatOut(tempMaxPos0));
    tempMinPosPar = find(AbhhHat==AbhhHatOut(tempMinPos0));
    
    if AllDir(tempMaxPosPar) ==6 || AllDir(tempMaxPosPar) ==5 break;end
    if AllDir(tempMinPosPar) ==6 || AllDir(tempMinPosPar) ==5 continue;end
    if ( forBack( AllDir(tempMaxPosPar),AllCou(tempMaxPosPar)+2 ) == -1 || ...
       forBack( AllDir(tempMinPosPar),AllCou(tempMinPosPar)+2 ) == 1 )
       break;
    else
       forBack( AllDir(tempMaxPosPar),AllCou(tempMaxPosPar)+2 ) = 1;
       forBack( AllDir(tempMinPosPar),AllCou(tempMinPosPar)+2 ) = -1;
    end
end
yytemp = yytemp - W*AbhhHat0;
ppos = union(ppos,AllPos);
end
AbhhHat = W(:,ppos) \ yy;
AbhhHat0 = zeros(N1*N2,1);
AbhhHat0(ppos) = AbhhHat;