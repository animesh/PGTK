#!/usr/bin/env python3
"""Comprehensive post-run validation for a completed PGTK workflow.

This driver validates source contracts, trace and failure state, required published
artifacts, all published findings, and runs its embedded parallel deep-audit engine inside
the exact configured Pysam Apptainer image. The heavy VCF, BAM, and display
identity stages use PGTK_AUDIT_WORKERS worker processes.
"""

from __future__ import annotations

import argparse
import base64
import zlib
import csv
import gzip
import hashlib
import json
import os
import re
import shlex
import subprocess
import sys
import tarfile
import time
from collections import Counter
from pathlib import Path
from typing import Iterable


EMBEDDED_DEEP_AUDIT_B85 = """c-qxHYjfL1lHd6&u<%udTtl$r&F<~-LRT>oZF6HwI+AiWr)9wq2}+1af&oC;icbIi`Y~?+r0nco<x&+3#LV=}^z`&=dK%CE@M2r!FRs${i)4N0ZSKpPZ2hN}@B1f7vI)FAiRTB|dU@}~+j&}g*}{vxS+?3NlQNlm=WixIK%HW{EDNv9lO$}lM#(SRG>1X}nC97fm8{Ed%cHu{XfVEbHyL~M-?4%tLCm4ytjzNJP&OE!U5t(gT%WruvTZ&~yxC1MyUn(xw=DobnR(aAI?3ZQ@yeS-0)*k5ZOSxT$IIH`^e;|^lj!VXg3!;l<p$vCF?c7*A_g#BFMNU359WZA144uZ@VB$?qv44d-ruq-Z+H@V{cV{+-4ejS&=z^N@;2!vS*GhGD)V@j9Q-jXi#zYrO|tg3>-3jxqUQU)|MzI{VQ@+~dLMsISKF1hNZ0dpeO-8qEQj$Gn`JB(B=jzd#5<&dc)2u&d>Go3{PAq`!(cS_e$4XQB=<IXHcN^kDU2onAe+TF=Hv5=!UK>;3dp2CdNY{7T*tg{tN5TuHgQhX#rR<I1pxUJ-xprGo-McYWUd!-G8!IF1g6S-o8X!ypMj><+4|Qc&%A8y4MwB0k(#bCv_-ruh`6?ESl_hpidmj+O6aCY=e*>}=Xh2EqStFA%=;o<0imwT7_go*kzC!=P{P&(stW+?<PM5g**smONnz5Cpl8?$tro3e6fL&pHcz6+OIMpLFTHrZ&dM0+t7x^vyIFCk9<%jqo994qA-61K$qG71#S*|w<?HKT(~W$-iHn<MdZm7rtEGJYxyaV?Ia4t6L_QTa+cH(P3mjolrZZKuxL1#5oG;R)YFVbMgeN=8mP?>toSK~NakgETNgjA}UX^(|E4lI}E^(}~;yL`}g7SWYbRpjN*Y|-p1Q78RblEENd)$bgEZ*RnwUW=X1lDChf1rCWE?{YATr-bTAouZoQ6{Sa&`x{d6UlVBu?K{aOZJnl1z99BKnILMQA5{z0Emx%2k1%gqaF>qofV;A6XA7IwtbUl`GM#Aolcm1-o)!U5LDjogs`UB@-At2T7psS0>oAJz(HoT3y5M_2*=?{JHQ%zr3a`f>w)W92>lS?wItu@d8_zy`zUB`Z1EXr?^YVWyBNKXfVGUH^T7x>%J>KD<p_F#ehdcF^B!om@8TtJ#-fCc3Nyer$mZxRu5RrhG#JKn2+R`b#1DLikS-jsk1fHOdi4Og6oC;B0*K4(6CxLf@BpKS5@cui_dY57N8kht;Ay2N=q5~y2xLmTvyX%=0ST0hFk&r^fNDh{Qm*KON$ShG5R_Wm1DA;oh!lHr9w-mg=+;;m0{;i*Ry6c#b*d%pYFm`vRpNo{R=nF}xwtL4bD(8lT!ow#v$Zf^VXZcl(`E@*w>Wh%l->lE0sxGG4km*WSX5+bA^rnvJsgiO2IH=`1hALLzAk|Z0_026DOA0jw%!kLxuJ?(iX;gO+XBu9cq_I^dkJ!5+4ZQ#ESo3daRm%LzLvj`FVZ8l_kCoXz<`!T7h%$?j(70NJI&S!&ydF##+yyDp0`^<Sc=^s|MRJ9A7+7n1G@}|0!8@8DG?)sLh<C3<5-~(|2QRL)KDm%a1*5<w?MM2XGvSU4RE3zFM&-$5-OZ;SDV5p1s)&WCih~EbSIDW)Dh_Yl&2-ANyS~8v4@-x^w8`RpiRmSckKdD8It9^5GxA>QwOu3((;Dp7G{9kw*LuMBUuBQKrr<DZMiu3E4*1JpHNBb;lhB}xdBGH6aWYa4h`25WPEr6L->(iw7H%b#5xB2_INzF;>&b}m3sa^O6>D9D1b#wg0<(T*Xt}#1ZUDcg!JbIRBz-nst}nU-~l%up~G4wvn7<a9{|%`tVgU2eh@b#81MsO@DLv*-WO?Jlo9d9HZ*g0Kmdwe8iT%QtODX$X#);X!w()HaZB7kP#edg_-9veaxtVpGj<#`Ok)K*|534h>K;qwIMnpHU5|UB=~5Nhk4B*6{lxbGQXA=m@M7aKPe55)bKRElLd;YkEG{Axg;cA{x{X2w6u^KE2@@&X8>sZ&>wW}Cgs@m<aoI-n5G<71py0#C0-?pa0#Cu!j4HtF@u3zCBb+Cn+gI_5Hxdm9AH4rDHV2}21~i<`(R@P0kotimf|eccx#x#haq4#h4}OQZm!e|3SfroZs-W6pRzQJ)S-`&P1O*oAoghWwmGK^xRnS|25{vZL#Czo()-H`|U@i|%NX$J?mInKvS;N?xmjMS(+yNSbwtNQ+C~FI0XJh{ZF3}**v%C?A+)=aiI?f|t(rcIv%nSAvY-ypYhW#b%{$xH*uF)a|S&<coL^$?n!TdO0#$Y8Z?!iC;Q}5s^F34hDBp~t75YDn3*cdKdd>tcKriwAbN6=7p1~o!7WW+k)IB*HJ>lMnHIf{0Ar=a41K_v$rL6eCwRbRIB#R4oM7(bAF9TY`lz<P7CmI9Dn>uaDw-0awbP*rX>%Y@9%%c(^+svXhhVc`9H;QbZ&hxg2TA9%gM`@2OWG(jHFcRj+Z_nmt4Ty$-dl{vMT7bZeH`wsC57(D?*l>(PHhY3p4OB&wEKXYIpOA-YC^@I8o4eIotIjG_`MLLM`<oY)a?TCgp1TnLo-DGIL#fy^s6v<thZ3~0oabesk(71a&)#U*IY#cM}ktD~CBks&^8W)n=w}JQP!263MBG8sh*S2`@g+?1G3R-5f+rpRfLT*S@f_AV1+W)i{E=j0_VK}YfTsNuQZXofZGQMv6Z-@RqNkOfMQ5pkp8`~S;{jyE+`v~=zz%J1HI-NHb(Ux~g7t8qCc||)OZLqaba+3zOE=ayASkP>#Y#Si}Dz26AoZTtsP<cQ$f0pZsf$bVGw4x_0RNd`r6XAM-Zq{=gGlW}F*7kvi-TEDSh-7z1ce8~6-;VR^?TV}-v52HXDV3Dpxw=)IK_V6%B}#<-2hvKSU$Bg-Csj{64pf_tz!_UtqSOQ<E3njq8nu#8wCKp%C{GaF1_V>2e;!gX>#Y$GHnMX=C1(qZ=o959>543o(4s1KnPnB~XMLV_&xbfq<8}EO)d?knNq>h`v3IRKG*lU($2L?8xg-=N9jg`T-84yL2o<0vN`TWHRSJcvF2pdXJ*cRoAq3hFfcDS(<FU`8uv(+q7caJfc-84n*O`fTL=zMknkE@1zUP~-*5Z{G#k5Y+JbjRcDBJm<tg`Xxe>i>}jRyUbNWNF5C-ol)Q&SRHpbbple-vICR9hvzszoKWvb;*QD~A<%sH+u9Lvv^bLb|E$38YQWg+q!dzpE;xmHv}?_iNBJJUtowJ$esz|M2vUi76k?Msjq}&poIp=CVqOEBQ}BAzF$(!*w7HeAUMDukOpF5H`OvNhlTMMNQ}_Fkxx>BLMm*v=qAt3Sgw+&sn-|+dI!N{(YL6PhOHXqq~;#kk3J^6}cGDpW+iiL1QMYmi8tHF#$*j#Cll#8jO^G6~9K%Phi(jO_HW;@pa4^EwW}_E%$&D`Tvh>JvnyU*pvTVff{j3Z+_XPvs*e2M3M5kKYZuwxkspNptl0bUZVq-Wl2`XOB7&vx{7mjTNFUO>$!N|ZeY=It@!sGXduOdF(u>Jv|Hog=xMH=QgW0p0nIDoDq_n<Me>Vk%}0LZ3j?}O%1cMvo6>kn;rda%Qv(I$x64quqpDV(cIS?^YxJFx!Q1x5WVo3j7$ovpAnWx|L}K}1sAzv@8Vv9eq_Y3?1T&$R227g(NSxLw3Vo@7D^&r9hMexo59GNKdra6xS?2>JHG<xeZpbvIkJ#FsDhZn;h;YArY&5s1_@ze_+&Hd<FynDsi!m;ea(2^3ye^|8TkyQ9ItA~CTJ1$|jOToCB%KGUpi^yQZ6LRuSy*kWH=EN=6;#{EZ7f=_w&jeHn2%MmdmPeD!qyOM-z^!RwR=Sw62@k)MJBpldo`9|ZTG@YYjTU)UL6Ke+cApF`ql>cFrLrJMO&%(3mQSimRZ_irTD9zJC;c;c&b=L+i}s^(PUoR<mIj=X>F6E`X+TAsH4F9R{Ei88#!33x<D${f#Z@xay1B0;*<~VPRYXuHH<F)(5u#OAX6ouELwq)s1ipKKTV_kzW+Wt9=_?1qW`)WjDCvV4Nl)o-kNId!hcX*JR0MmRxfDe2V$mp#PK9WR}9_W1y39QkNwf9J^Yi4^LNAJ{$vpK-wofKz8{=Uq8|o7jh#u@ojr=!6V1A`^A*dVG>zYh@rFqFjmBGp<xw>}P{DM}j*4CIRI@-C1rCG=wo8;P!T5mw(2vBm>x{b<Dj(6h@Lu&iduX!6v|F^&;dk|OAiI?^-|yPL1L0sMOTh1v;57{^mSlccJfJltB^V$>KH_Urd_j@<(*USw5;{}Q(46_;X_RfCK(Jo*`ow%~xcWr#&@4`&l_B>?O8Q;njnQ0f_{d~ZP+Ptai~&cMm^6n^zYVa$3VgZ|D8gbQ5lXk$6cqm$FG)=|3`AW?>)w?^YBVZcBatm)Q1M`1*BP#?3Se?gMxp+rWgzP)lM$yg?>2YN0Sq;ELz7p{R!H@D|Bn<>sePF)XYV`LLag=1J0}pLIl8l9)`seCm$fV1bsDrcNi^vku#=YeN{I`c?sh44&Dxw?%uA>;wfLy}pRXPD<fO4&_vBghgXYHeSP6jn;z=#Z*Y9_=n(h>Jc0e(-LSO8+cAdJ|?(Bjhj#w-JxVZejvF&5+Thf8^WR=||^=h!8uQe~A?VR{migkJJT*E?`h608z&!e-|2Go?O>0H`%ay6Vt+lRfE6>^yFKkvP?R{u$dy}!p3-rw8d?(gzkmsG92$FsYSX7?YkKJ$)3Z_GCXBx4pPl#r(Af-G8mc_2dLg1us^c%3ei0v)s>F7gFvVKQ1@wk8K=5tg4zzh%g0Xir8HgK*Cv%!&w2hMdygz!ym2sxj^hlCf<Ovx+tFhzNSHT1!8McnRC6K)K!apZT3>OB?zerzrJ%&*$ky@B632*MspypE~TEFIt@ORRI|FNr7q00e3pB7EXkO4ee`g+BHuj6-u#|;D;E}NoR-Ujm8)Uc%mD)S9ctNNDwVa&u15ev6$Mf5yu=P!$38pMtDL%55}Pff_=0Ol0*kECmu-8R5u#1`Fi+nU|@3LtJ^Q)nhH)&`h^g5vf7jwXp7@}6e#5r%a)FO#4uSPLoq*z@rZAHKm}~4=N@>MJZn8)HP2}W_r}Z68<E8X3x#YSghWd{8Jl4duhQjxG{ZoDP=ZCkzdXDZKVh+O-D!nk-7%~)nag0QFuz`AS8e}~>Nr3k*IsUoS0|7{9c{=6E)n-_4(w7KjBwo!Mj^w?Lr`*+JHwvbzJOZyWXvwzbxP)h*a|bTT3MDj4ftRM&&h;mo14*F>y%7}$W};19uKJr3-YiLU{VNACoo1EsDVt)9;i#Jo5)%Q?ulw$KKNV&Q{hlk;vc>U@qDa_@f2iZLli-uov8pTp(=n-(2+6Pmugg8<dHOL+cPKf|5Hd;vygGp^(#IEY|(S1PM65GvWb~#x$f>39)rPzddqlqHIKawtazQjp!IWhg_cXh-2%QnCNC<83kk4cxRnh0!4c~y{1?DvoDn^6{YM&@ApvAlccPaEM^oIF6$B7(rQM|%Bh~2|aWr7Q(e@0BB+kK<%iYnm!{?oU;86YzCPuG9+Zs|hAz<Nko^3bnqYnEY5j?wcw(Hw<_6bxEu~5Z!g&vU|Y>De}sCe@^n4FZR2bV9A8O$XVE$#-2*V?<fa~W46qr^uX$StG9?U3+}zStsSyNy{W;7D3q1nfj<<HqJP+`INr1imi0d|L7K3Ya1m`oTXwdw>3JFc}a%9`}#m!qa1oAJq_olTjaJ;P?HLXgujpF2*DQCxl#a*{YD;Jn-K~I%!+BX^0()>$kkPBB;7%nsHg8>;x9z9d+=3@0F`Wis7+WSB`E3+~Ozt85g&nfC4HF7<t}N(`Xcn(TkxNjDFSGp>1*JhK_ZbsY{@O{uV&EbD;w0SmSyIES8jTSj&F}6-fR|Ri40qEmEoSUjbDCB>oE=_1n;U9jD7}o&ZYb*Gb;uKNg<??oCyqWp-Wo@Iq{SvxHhP-$Ff?f7WeBQS|iv$PK=FL#oj@JJDTd3^%wK4Wf60lQ)A=bb59=cz-_miFSgGBD~S*Js<@cfIE-kpo8d;3pzD)C288}&Y&-la0kUEnWc+#<}K`=O&#9q2bS{BU^j78&P|e(#1xay8@{?M!gX4Ou3*c{OKmfWI^m4;4=F+UT9#WrZ8b#)R>D*2yH=`YO_Yh_&eS#PXy>3vul+{<C3n7`hKyJQ9b1UAs$*na1Ne2v<oae;>B8I9ZrZZlOhUt!D6J!#aUXg6`0ZdYi7rn2AApMb-@l`)s31O<Y*6Bj^?Zf=$t3{aP`czR3Q`gD0c7--kmA%}&f6Kr<8rJIF{(gbHb21o;>%^S#GBZ3?HJ=`Q8BI2y`&)#VtD5DuedL<wl-Gq7n#3gVVl?d@ScwZY_Je3nwb4?946d#uUG}LE}t6U`pCNBD&Dj$B5@;tZ5{toM;<zy-ZQNzHQuRi?l`^V?4>nhiB9A8VDh=dVrvE(cOjy=ZK$qv*9-m;7Q(c~$iTLy?I~f^niVoCFmR#-jkrL0NTwtBrUrolbS-RUS9Up%M#Oe-<s5}++n@sPz%wW$zz_l;(|-@WH$YeM()+&u9(|k1buL{+bO_DMxr1c!C}bbWLLH=6c~SsB3*laMLJT>UPMpplJz!WMkb~-kNx;y$=@FD}*0S5|#T6a>4;@2C@lwCh&ROvW8Ub)r&GCp+E9Of1Mojo314D$DFA@+>VoYpmYqB&ac-gkq`M^@|UxHHM(XncaFsQdN{%YLD4j%sX+gJyNe$_Ts<89u?a`s=hjaBO}*~X?Ki2%&KECEW5ud<)jPIZsd1M3XDdhjhUU><h|`heRy*sdIB2iC#0b6hO}n`hO5eKuWpCLJ`MKTCj(N6(s1S^ftkd&Q|FBJJ`gaw=mJNKW<F(0iA~b8oY~0>$G-XIUv)xC>KqlGiye+J*+6Uf=noZ;Phpj|imvGRt8l$>l9!g0L5U56mZn7rB>0{7N8ujO)75NvE(VAt!IO*5%<;Z*$Cr!YZ4piqfU5xwu`czey*1wBD|+Fl&kQxF_&}5(Oy(;K<pc)|t3Ep`fT1y4b#ibt&2iul9)7xEfDP2oco~C8jdnG%<s{$A^t}$?)_9EQ9E8<Fiv4%L69KJ4C{RftN>_3?&dQSApx!F5^}JA4li?(PTIoLKlW5V$zLadKeo+k)zx^@&!@;e}Hi*5(>q@!){hd2_{XtzT=BN1&JsiytqVLhfbO#ZC0H`HGv#c+iL_<F6hK#KXn`Wx~kqX){)&>4SN3oMKyc-SZb$p!X-E50u2X$iO@{`0sl^goxyrkb5Cr~X^l%YoZaNvDl4+pj>gAWhxn=WK&{ALxXBf~O;*VL-0Zj+)2T2HMn0TpEDF>g0-{GL8A3F}o6kFjDl5_rX&S2ClwQdPfoQI)P$?a24PBs__W`8AoXBjf+#z&Ycc%&>?I=<1xQ@(igrR|r9h~P0l3Fq!QbiV#;#HKYcj8Gx1r={6&IhOn_LpR(viMS(VJ=^p#kydsF8M8X_HevR;`K;smkQ7fUmKw{s<&!YJPWq_1U^JL!}d&#D;ZPc7Y?u)up!h#*ju)HE~0AR(o_@xTY-_4M^<lQtM2dXf-2fQiN<z6=d36hHL)-4_poO7haQo{b6mC^zA7OTj{EP({EvpGANr$V|8!zbO!tUjElPq8>epTKEyIZl!?e(XY}<PAhn_o+=ZXgQ&*I{A_+J;!99-;T7d7*G?zvRdGy#{2tTW+7&WT*={qlgc*P=pyl_BE+%*o9rn?^AfMKf7uG%`K<<gN$)A&S&tk;+A-y!Z%Y9rM2gah<J+Z?Ct@CH@WfIraa0FglA)2B&B5ho}9?*(iEF>K~(wH;KfXgFZsx=#;weYFm2!cN2{fj8jK7E1lIL*OHm?oTBY8??5o;K=9CxafIQ>?0pBkAUnI}D4bJR7Lr{VNabWyA~D!CWT%@rDWd~lA;sL*ByGB`Omjydo%dqK&ghx<cLg4tQlKUuQ#ZQrffI00pYvXBFp7!bM|aF(jcM2>`8)M^FGDSJQ@f*xKlx;r?mEjh?Ndb+sNcor8HZlnikPQU53n8r&pHRWJemdu`tZ%mX<%4G@_fWH8tTufZDLg?CDc+Km(-{qS~1YVT&ut0IOXUt^xk0lQqXgy7nbgmsmz*1BC;_;ZJngAcgSuFunF$Ko~_!7q1P9bd6$SeF(go;0g~O9RH@PLy9J%S^gi8S^h%Bq{%eX~QMpt<qotgd%X^CZNecAd3{Q%%#o#d(d8d(Rq$B(ecq&H57_VNjRp9J$nFdbuAU*PKEyHe+RbkJIo0j6FLEo-qI;-{085TpsjbDrg*10g-Sf_HF>zWtFnJ)Hx*(RrrRMZ3Z;V}8gUs9zMKb0;xt(c7)ykm`Uw=12n(Ln5r7xE?Imq*Ysj-W+!vLJd<ZUpPpxEsVYFk0CPi6Dm*PmWor96s_(H5F~ts6)5476e6yVi)zKY{|7rXSaoULAbC35iM-=+|y#&`qFSfxV;0uUQVa2=38`iK-j^PdVloVF1!Ogw|?QZ&5_|dkZp+mlmhAE^>t#_Fb;f7y275TA)gYfQX??80BK;f@Zu$)*PxDW&D~F~*XW*bHyPa^)GJGg9Cb!carfb#&EyAFj|9YyJ@Z)E%vN*j7)?~tdE>z7q>H+9gRu3G4B}=74~>pObW}r&KB*t)5}xj&2I<a;;y#R4uxShfCB>1J6q^iej052n5jUEEp97meAP0lW6s-bci(~~=F%$4Jap-$M*TC>F1gIki=}=~Qe4TjcT`xk9d=}m+K~GU}lWuV8!12;GXHzb35`d&HyQyGS@7x+2sD>|e{w9n$xaIE1cZ`=EOP#_QAnKH6p*IwB@_0^ljYZZonU-ACl*j;~#Qqk2=D^1qm@?H4SlH<<m1iO9SjPw*HUL5kSwVY|6=Fc2L2qh{!>akny3$?@m#c7ago**>=%yxVg_n$klyij>yjNLc-5I!t=dMX`G8~`3>;F_6JG^(YmFNSvcgxl2_<r5GO)lWQxjS<@uQzgc?hXPm+KGthS+?A+*2H;CSDyB@9gFK1Y{m&>{E*!HTD_EK5*k(o{DHXNMd^riLQ!!O#=5I(1KVxXRaS6({C4oZPd?%IB05<D@*V`6bb9DKHZNeagfz33(yXdtM!G|VmrN!tHPP&MY6Oqz233vpg4DH<cVygF=Hk`0<Evj*6Kjndnqk<woE&+rdB~gxMT`W0V6xVRtJ=oCzidlj)q8`e49eU4io1X?qZUnv&!p5EZ&OS}3I8;yw(73BW?uVP2dNg(Xz&IP37ly8Dx<6L3bU&!pRfkPiGB0mDRYp{bma{ntmnTA&R{LjZ{Z5gQvD{5Kt}~f+L5YvH=(a)r-3q6S0Axlzb{1>-}S!uQYmg~?Cm=NqrSNq5Lc<WW@YbCjxT87e8CP5r;WhcQP^M=6$+q)fdfotfupUWslyis8*jWe<+8P}l6E;UxOs<A9?D1fl?%F-OCP%ki<|i6pa0Ts*JWOz3$2}yHa8G(T+Gtc7dIRWolCnD-Xx#r=`}FTcIUGD^4qEYmH_UncKt<J_o}dX+kCTC7h!hGPEbi}DM46swJvTtv-bb-`Dh_7f?#+Gv1A^1N6%Xvd0W$z_d8*At|R0d&<Yj;%+EzBQqkw>^-ZZ)0#hMa2{o~*xcj<+xnjC6Fx!D$;;lVwTUW^aX27uXyQGNX^?f9@?k}Z}KyHcd4yoh}ov;6?i>Tr@puT@feg{G4heYFt%yhV;mUE@Iia+N{Gy~j~30AX*L*y*nvt6IcWRdsWj>!iMwrLt7|MGWS9aRK=4L`kCvM~jW>Dqx&2V;>49rdZrRT?Vsn=_-eToYTEK8Bes{gIZbtwa@PG^^J&cG@KfdY5t&aPkA=z9dIM2{)_?2Uf~fqY-Y66}WCyjoNO|S%aEX(P_z>)O@(O)9|XzxFA_)eT17VE#u<0&>;rx9rG7L*qVm&Mt6<Wf@W+UasMh*b;@<cpU$xpTQ!*wm?6fz@hBCA2ZsL`KFAIpDG*XFXQ%$opM?V{%Zv<Y)!XfS)n6cjmI+;$>B0u}<q?B>@6lHR2<!8R1ol9xch?g;ef&Ozb226QnzNh$-X*@f8sC16!BeNj$oOuVQE>03Xw_sngLXC?sh~2>I_|7Rv!~B{cmkT~aPpHyZFGS`fa*Ldb|8IR_zf>MAjEgqV9jCM0nPeP!&%MOD8z(pW>c54gPBAXf3$lk{*^x|q2iz<oe+R62a@cG6m%F0M+-J^W*%t`TV=Vp$%-R7x;B`6pnvKk3H3@j^3@x1=9~9?xmAUW$cqSUZSZx3=;HKZEUzKRh_5BN_X?@v%E31<556uU3~9_4;adUGcI1okvJw_;egZ^n3R}_;2i+;5=uD|&UD-=^U82_6E2S#>l-hJ@%t#$tInXk!3)<PXol!pz82#n6<DP=n=c$pzPQXpNuNN{*`<-fP3IU}tpJFZJqGUq4<JQhxKW=&w%5yVycDux8eeOeJZ7vy{k#lD{t>#*!kP+!(>Uj=1cAK;RIV!1R4CAqqFIDSkm6S2)C~+k7K~I~KKBl|0b5{>nDtPCK;}D?d&HXneoveVOu<{+?Xey+nBMU;)inv$3^)8{)0ZF6}6L#TCckoqKfZ6S*oI-jCPJLH(Xc4V<C+QJP8;=Yfz*;DyfE4mrlee9K_o)=hD<V?Y?pXd-HOZBX(%nnC`}F~q8N<`bVDzE?E&}P<Kj}~UmRfL{sq-FB_Cb8?L81ec1)A!Cb1Ih71%s-bF>Xokpwy?1Ug<Py7Q^xT{^a<rP)oR%Z)$t_X8P#i7l6=&R$m&JXGE={4jPWmu~+S+s-#!D#8!eNAT(8=vJ{ywD>-|nO-OW-3QeD`3yv`@=K%inPov|r(@CF#(>N#%@TasnzGN%WY4?YCl<50(lNE`$@pbB9J^uD9u?*w<v5_>Gw*VE(yIFFM4LOQTgxL(D`Y9U|^3zdteg?$l(3OK18??rTrX*w`UxWvBRTNX5t3dus3QNPQ`iF;|D$A*;WSsKB0)X*r0N~|eWlI_)`<?2)PJx&sP5Xnt_m3x0|DE217pEtK(ec^E>BN8X*!<Dp^<XqO1tP2gvu}~;dyKFkXL{(B9sQu^^=M2ChDZ5;IRN!RzJA=rv<wWr`YH^aefVb}NJIIj<LA}a;YSAosDx2=yrOXdqia}Bl^9VS*{O{q8MCEh5Vb8OBe)tOJGBiaV!Hw-@H8;uOphoyU7S~TwWjUHe{RS_5^;d2yD@b+2a3A;&R6D-At$>^x#(1S6Q)3=wxNR<xALrN;w2R8C63~fwI3+gZ!58{hMej)h_Ap%uGg)RWbT@{(yY5f&;^x5dEG`~=kz1xgvwyKPk_)3OIP*k?o!JtGpkc!y8Bcd_iM27lBB#uN!Ag9*{<40Y?0Ge1nbeL5b{9==d<GquP~8|4hycicrKkCyb;Z7t9l~`A~M*yJEc4b{?1Tbb;(IwMw%IM0M-uTf!e#KVV4B;+o7KN8;D+XXHl|A$UzOER*5mFkXbDeiO%ZRKWJR?pk-d$x7G^3yTCpub;VD`95yte>`}MWL29J)(hJ?e@VA$t_jdCB-3wla7t3^g%W+^u>!w^S&0Kc-bJT^{L?x;dpXm`FKa<HAtFZ*~z)>V+75+o{nOkfmZ)LPt!KUvQMl#PU9`TgfZ9@4?xfO*e%2#d9^>YQt@Y<`apbv+kVvn)-9?){xQLBid*tASJk4jvDliw4-BnBI#5sVVR%Qzx*&&Y9Ojm}J-gHeH~8jyhYR8BIo>ck@jY67yz-NHHmx%pbIet|2-TKn*H9BGiE7u^j5c{?0W&VWL|N~&w=xEA><*g}J#u5|9Hc&#GUg$EyXiWu{gsQa#zdXc{BBI)Pp=>kPUfUnILc!SXj*LRrzk0K8xLDW=JG=Y`(6UZ)5KjTLx1SLfB16ofan;NK%i3tgtCdXFjt<4%DQ4zS)4H^f{sV$HG6nd}w<4Iqqo@p%r2k{rOf<+9w0&TR$a)7z73SX4p2Zd|#X2YN*;u|XDkxW3H29mZ1w&2ZyHLM7V!u&o}#GhF!=xcU#J2f)vncrt@avzJYYU;3y9p7K%!P6#aj?1jD%U#Lh<sJRHa^g!^-m~)IHn~)>&X))FRFdkx@*BLNsK&Us(w`U#7u)4h=4!d@AN+TG@arGLgXwdMz7e0Xa+-H`d`hOeoV6a;-7Dz+cO%2zaMqVc-6<y?Zg*rNhUPB~NkS8!(Giev<+z};qeJ$!^K*T@J+o}~E4y=-tlZAk@ULIG9m};hqeC^r9B?=5n=C3<d_{&?-@%a^z3D8D0qZZ`yoD|XBTnM@PyzF3eBOo`t|+KVQq6xWX_>Xq#QZM;IcPk{XX^_fvWISUR-^(sNTx3lIPdDr3Lc*+ob-rI>7)VIj;eMDLMnYrIV?auV@0)Sk(wehHyQ_fKt%cdrT7>RhuF|7+)A7XsMA8#rm-b`wx^+?QnWdH5b4w$&1($Ctrj6fqps1rQI%Cds>P4Gh%?@Zk1bOT@xf(XC0WKByx_C%up6jsG#Q-GA-8HwA!QNhyahX=<!1*Y`7VmAY;`4GWwGjCC!m#KIwx!gs~AOYwJ73y`J90kmb5}c;&Vz&K8$cA9bQPvmpSK@L*nyEfv3nyelZD1MP>RH*Ge&Mir)!EZ)kr|!kPJYwJF+SQUQIqXkGSR2Bc?@Z$*5t@VUpIh^S=9X|(|U?Hp!<;(+bFXC5XY_U1{l;cKMQ5EOqgnK@;b?s^{|gevp+Q5ABweJtZ7kVX+^TOOWZk;qJ}i_Q_}SR$9`8*4qknVLO_Rq;x&(zn+)UbV5ia8~;ti3VjVw}j>)Xi?RYll7lH^UeSW>d?|r;9!HeZx>IL^zNZU?{BE&o4`BlJVFE2J@DDu>R3>zJv@DVR$ChjxOUBu$1_mE8rF?Jc>pj4wiuk5G`lhg&fUv*jO=DQ>usgchvV0k>X^i@QW0;NxfL-yty1y4lJm-KCf`?@ia3%&UlT^wD9`+e%w1Rwo{I|T69I%A6`pWZ>;z`}G#VFih`dB3%QTAVqO`$&yE^2Bj9uksGwOtB5W@3nm95QCCAs^hKwkI<uo<x7<42eRHO5dVty;M)=!li`=;5_-y2J<#R5g*W&bg~^eUJ^Ox<YXn%X7*N-YsW`R(+Gn4c-XBRr5r)9q6Af%_(y)@XGjVyNq)qMO1(hUoS8rDW{`LL9M+H#CacnF)=M*v@>`mKC`%#ly*>N@`U&RQ$+Gd0ieYV9uaQoj2(G|BtX0dDZ*JpZ&-TsEGfj-Brs4g@vd;s02lgOd>4;rV6s5H{`%fCJ0XX~3J=;)A(ZF^1uLfyj*=xfOX>P!1|#|zioM~R4??-1IpW7cim(Zd<0SDuYLWFZ@IE&9R<X#+5%`gD^wCNJ`jK}P3o_a3w50{wXmEZunuM#lEzM9m$pCkzo1HT2nRiUSz3r+H`D8hp5+<Hc?~>L`lt$unBcgNpx)71|ev2~!o3+Jtp-)dqCs*S?fS#(Kx{wFB)~7D``#0Dc)7_XuuuTSz<*(j_(Ednrmb9J7|C;X-6TaK1!&4<-VGja`{(vHkzw`JHQf&k>8GSloGqlIw^qrCr3}U*riXrM;J>)NqNFF2KZrA9Ng>i8X{W#ZGvqTk8fIf8rRzRzXpSUNH*ZYqO+QoFxScH)T*awcLCx_LMLvOf6BFEH>JXOV<wS31j7ZhYsv93ftBZ93)P2{}EDJhiis9oyfGSrb(l;))rFQf4FlZ2wtKuRu_a5_HpO9Y_-?n#P}QU3{wkTrE=YEYkd^EZl>d}G$?Gs3$H^=9GtZ66j*m+><BoR46PGoFZC2+oL<)qN<x;Y1%MAK!zxSPed>W&35P^}o`fRO$"""
EMBEDDED_DEEP_AUDIT_SHA256 = "fa41337662057edb6dd20b5b9d0a735f3fa26bfe1a3179efd85e9a9d5fb09539"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-dir", type=Path, required=True)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--job-id", default="", help="Pipeline job ID; auto-detected from the newest results/pipeline_trace-*.tsv when omitted")
    parser.add_argument("--workers", type=int, default=32)
    parser.add_argument("--pysam-image", type=Path, required=True)
    parser.add_argument("--nextflow", type=Path, required=True)
    parser.add_argument("--host-python", type=Path, required=True)
    parser.add_argument("--apptainer", type=Path, required=True)
    parser.add_argument("--pipeline-mode", action="store_true")
    parser.add_argument("--output-root", type=Path, default=Path.cwd())
    parser.add_argument("--max-events", type=int, default=0, help="0 audits every finding")
    parser.add_argument("--reuse-deep-audit", action="store_true")
    return parser.parse_args()


def run(command: list[str], *, cwd: Path, log: Path, env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    print("RUN  " + shlex.join(command), flush=True)
    result = subprocess.run(command, cwd=cwd, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    log.write_text(result.stdout or "", encoding="utf-8")
    if result.stdout:
        print(result.stdout, end="" if result.stdout.endswith("\n") else "\n", flush=True)
    return result


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def process_name(value: str) -> str:
    process = value.split(" (", 1)[0].strip()
    return process.rsplit(":", 1)[-1].strip()


def index_for(path: Path) -> Path | None:
    for suffix in (".tbi", ".csi", ".bai"):
        candidate = Path(str(path) + suffix)
        if candidate.is_file() and candidate.stat().st_size > 0:
            return candidate
    candidate = path.with_suffix(".bai")
    return candidate if candidate.is_file() and candidate.stat().st_size > 0 else None


def main() -> int:
    args = parse_args()
    project = args.project_dir.resolve()
    results = args.results_dir.resolve()
    image = args.pysam_image.resolve()
    workers = args.workers
    job_id = args.job_id.strip()
    if not job_id:
        candidates = sorted(results.glob("pipeline_trace-*.tsv"), key=lambda path: path.stat().st_mtime_ns)
        if not candidates:
            raise SystemExit(f"ERROR: no pipeline_trace-*.tsv found under {results}")
        match = re.fullmatch(r"pipeline_trace-(.+)\.tsv", candidates[-1].name)
        if match is None:
            raise SystemExit(f"ERROR: cannot derive job ID from {candidates[-1].name}")
        job_id = match.group(1)
    args.job_id = job_id
    if workers < 1:
        raise SystemExit("ERROR: --workers must be at least 1")
    for path, label in ((project, "project"), (results, "results")):
        if not path.is_dir():
            raise SystemExit(f"ERROR: {label} directory missing: {path}")
    required_files = [(image, "Pysam image"), (args.host_python, "host Python")]
    if not args.pipeline_mode:
        required_files.extend([(args.nextflow, "Nextflow"), (args.apptainer, "Apptainer")])
    for path, label in required_files:
        if not path.is_file():
            raise SystemExit(f"ERROR: {label} missing: {path}")

    output_root = args.output_root.expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    output = output_root / f"PGTK-complete-validation-{args.job_id}"
    if output.exists():
        import shutil
        shutil.rmtree(output)
    output.mkdir(parents=True)
    started = time.time()
    checks: list[dict[str, object]] = []
    metrics: dict[str, object] = {"job_id": args.job_id, "workers_requested": workers, "validator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "embedded_deep_audit_sha256": EMBEDDED_DEEP_AUDIT_SHA256}

    def check(name: str, passed: bool, details: str, critical: bool = True) -> None:
        checks.append({"Check": name, "Status": "PASS" if passed else ("FAIL" if critical else "WARN"), "Critical": critical, "Details": details})

    if args.pipeline_mode:
        check("Pipeline execution mode", True, "final Nextflow validation process")
        try:
            import pysam
            check("Exact Pysam runtime", True, f"pysam={pysam.__version__}; samtools={getattr(pysam, '__samtools_version__', 'unknown')}")
        except Exception as exc:
            check("Exact Pysam runtime", False, str(exc))
        check("Parallel worker allocation", workers == 32, f"workers={workers}")
    else:
        preflight = run(["bash", "validate_pipeline_commands.sh", "--project-dir", str(project), "--nextflow", str(args.nextflow), "--python", str(args.host_python)], cwd=project, log=output / "source_preflight.log")
        check("Source preflight", preflight.returncode == 0 and "RESULT: PASSED" in preflight.stdout, f"exit={preflight.returncode}")
        for name in ("test_semantics.py", "test_container_bindings.py", "test_igv_event_identity.py", "test_program_interfaces.py"):
            path = project / name
            result = run([str(args.host_python), str(path)], cwd=project, log=output / f"{name}.log") if path.is_file() else None
            check(f"Executable source test: {name}", result is not None and result.returncode == 0, "missing" if result is None else f"exit={result.returncode}")
        worker_probe = run([str(args.apptainer), "exec", "--cleanenv", "--no-home", "--bind", f"{project}:{project}", str(image), "env", f"PGTK_AUDIT_WORKERS={workers}", "python3", "-c", "import os,pysam; print(os.environ['PGTK_AUDIT_WORKERS'],pysam.__version__,pysam.__samtools_version__)"], cwd=project, log=output / "exact_container_worker_probe.log")
        check("Exact-container worker propagation", worker_probe.returncode == 0 and re.search(rf"(^|\s){workers}(\s|$)", worker_probe.stdout or "") is not None, (worker_probe.stdout or "").strip())

    # Trace and process coverage.
    trace = results / f"pipeline_trace-{args.job_id}.tsv"
    if trace.is_file():
        trace_rows = read_tsv(trace)
        statuses = Counter((row.get("status") or "UNKNOWN").upper() for row in trace_rows)
        bad = [row for row in trace_rows if (row.get("status") or "").upper() not in {"COMPLETED", "CACHED"} and process_name(row.get("name", "")) != "VALIDATE_PUBLISHED_RESULTS"]
        metrics["trace_rows"] = len(trace_rows)
        metrics["trace_statuses"] = dict(statuses)
        check("Trace terminal statuses", not bad, json.dumps(dict(statuses), sort_keys=True))
    else:
        trace_rows = []
        check("Trace exists", False, str(trace))

    source = (project / "main.nf").read_text(encoding="utf-8")
    declared = sorted(set(re.findall(r"(?m)^process\s+([A-Za-z_][A-Za-z0-9_]*)\s*\{", source)))
    observed = sorted(set(process_name(row.get("name", "")) for row in trace_rows if row.get("name")))
    unknown = sorted(set(observed) - set(declared))
    not_run = sorted(set(declared) - set(observed))
    write_tsv(output / "process_coverage.tsv", [
        {"Process": name, "Declared": name in declared, "Observed": name in observed, "State": "RUN" if name in observed else "NOT_RUN_OPTIONAL_OR_UNREACHED"}
        for name in sorted(set(declared) | set(observed))
    ], ["Process", "Declared", "Observed", "State"])
    check("Trace process names declared in source", not unknown, f"unknown={unknown}")
    check("Declared process coverage recorded", True, f"declared={len(declared)}, observed={len(observed)}, not_run={len(not_run)}", critical=False)

    ledger = results / "failure_logs" / args.job_id / "failure_ledger.tsv"
    ledger_rows = read_tsv(ledger) if ledger.is_file() else []
    if args.pipeline_mode:
        check("Failure ledger deferred", True, "written after Nextflow exits", critical=False)
    else:
        check("Failure ledger empty", ledger.is_file() and not ledger_rows, f"path={ledger}, rows={len(ledger_rows)}")

    # Published findings validator.
    explorer_root = results / "igv" / "findings" / "finding_explorer"
    findings = explorer_root / "partitions" / "all.jsonl.gz"
    geometry = explorer_root / "event_geometry.json"
    config = explorer_root / "explorer_config.json"
    server = explorer_root / "server.py"
    validator = project / "validate_published_findings.py"
    check("Finding Explorer geometry", geometry.is_file() and geometry.stat().st_size > 0, str(geometry))
    check("Finding Explorer config", config.is_file() and config.stat().st_size > 0, str(config))
    check("Finding Explorer generated server", server.is_file() and server.stat().st_size > 0, str(server))
    if config.is_file() and geometry.is_file():
        try:
            config_data = json.loads(config.read_text(encoding="utf-8"))
            geometry_data = json.loads(geometry.read_text(encoding="utf-8"))
            check("Finding Explorer geometry count", int(config_data.get("findings", -1)) == len(geometry_data), f"config={config_data.get('findings')}, geometry={len(geometry_data)}")
            valid_regions = all(item.get("regions") and all(region.get("chrom") and int(region["start0"]) >= 0 and int(region["end0"]) > int(region["start0"]) for region in item["regions"]) for item in geometry_data.values())
            valid_structural_roles = all(
                ({region.get("role") for region in item["regions"]} == {"BREAKPOINT_1", "BREAKPOINT_2"} and len(item["regions"]) >= 2)
                if item.get("event_type") == "FUSION" else
                all(region.get("role") == "JUNCTION" for region in item["regions"])
                if item.get("event_type") == "SPLICE_JUNCTION" else True
                for item in geometry_data.values()
            )
            check("Finding Explorer structural geometry", valid_regions and valid_structural_roles, f"events={len(geometry_data)}; regions={valid_regions}; roles={valid_structural_roles}")
        except Exception as exc:
            check("Finding Explorer structured metadata", False, str(exc))
    if server.is_file():
        generated = run([str(args.host_python), "-m", "py_compile", str(server)], cwd=project, log=output / "finding_explorer_server_compile.log")
        check("Finding Explorer server compilation", generated.returncode == 0, (generated.stdout or "").strip() or "compiled")
    if validator.is_file() and findings.is_file():
        published = run([str(args.host_python), str(validator), str(findings)], cwd=project, log=output / "published_findings_validation.log")
        check("Published findings validator", published.returncode == 0, f"exit={published.returncode}")
    else:
        check("Published findings validator", False, f"validator={validator.is_file()}, findings={findings.is_file()}")

    # Intended published-result contract.
    samples: list[str] = []
    sample_rows: list[dict[str, str]] = []
    samplesheet = project / "samples.csv"
    if samplesheet.is_file():
        with samplesheet.open("r", encoding="utf-8-sig", newline="") as handle:
            sample_rows = list(csv.DictReader(handle))
        samples = [row.get("sample", "").strip() for row in sample_rows if row.get("sample", "").strip()]
    expected_progression_reports = 0
    subjects: dict[str, list[dict[str, str]]] = {}
    for row in sample_rows:
        sample = (row.get("sample") or "").strip()
        subject = (row.get("TK") or sample).strip()
        subjects.setdefault(subject, []).append(row)
    for members in subjects.values():
        baselines = [row for row in members if (row.get("baseline") or "false").strip().lower() == "true"]
        if len(baselines) == 1:
            expected_progression_reports += sum(
                (row.get("baseline") or "false").strip().lower() == "false"
                for row in members
            )
    metrics["samples"] = samples
    metrics["expected_progression_reports"] = expected_progression_reports
    contract_rows: list[dict[str, object]] = []

    def artifact_contract(label: str, paths: list[Path], minimum: int = 1, require_indexes: bool = False) -> None:
        existing = sorted({path.resolve() for path in paths if path.is_file() and path.stat().st_size > 0})
        indexed = sum(index_for(path) is not None for path in existing) if require_indexes else len(existing)
        passed = len(existing) >= minimum and (not require_indexes or indexed == len(existing))
        details = f"files={len(existing)}, required={minimum}" + (f", indexed={indexed}" if require_indexes else "")
        check(f"Published artifact: {label}", passed, details)
        contract_rows.append({"Contract": label, "Files": len(existing), "Required": minimum, "Indexes": indexed if require_indexes else "NA", "Status": "PASS" if passed else "FAIL"})

    raw_vcfs = [results / "vcf_raw" / f"{sample}.raw.vcf.gz" for sample in samples]
    artifact_contract("per-sample raw GenotypeGVCF VCFs", raw_vcfs, len(samples), require_indexes=True)
    for label, patterns, minimum in [
        ("per-sample non-subtracted protein FASTAs", ["combined_fasta/*.fa", "combined_fasta/*.fasta", "variant_fasta/*.fa", "variant_fasta/*.fasta"], len(samples)),
        ("progression subtraction reports", ["progression_vcf/*.subtraction.summary.tsv"], expected_progression_reports),
        ("IGV BED or BEDPE tracks", ["igv/**/*.bed", "igv/**/*.bedpe"], 1),
        ("IGV session XML", ["igv/**/*.xml"], 1),
        ("IGV batch files", ["igv/**/*.batch*", "igv/**/*.txt"], 1),
        ("IGV event/sample manifests", ["igv/**/*manifest*.tsv", "igv/**/*manifest*.tsv.gz"], 1),
        ("finding explorer", ["igv/findings/finding_explorer/index.html", "**/finding_explorer/index.html"], 1),
        ("final MultiQC report", ["multiqc/multiqc_report.html", "**/multiqc_report.html"], 1),
        ("results catalogue embedded in final MultiQC", ["multiqc/multiqc_report.html"], 1),
    ]:
        matches = [path for pattern in patterns for path in results.glob(pattern)]
        artifact_contract(label, matches, minimum)
    write_tsv(output / "published_artifact_contracts.tsv", contract_rows, ["Contract", "Files", "Required", "Indexes", "Status"])

    final_multiqc = results / "multiqc" / "multiqc_report.html"
    if final_multiqc.is_file():
        final_text = final_multiqc.read_text(encoding="utf-8", errors="replace")
        check("Results catalogue section embedded in final MultiQC", "Results catalogue" in final_text and "pgtk_results_catalogue" in final_text, str(final_multiqc))
    else:
        check("Results catalogue section embedded in final MultiQC", False, str(final_multiqc))

    current_html = [] if args.pipeline_mode else [
        results / f"pipeline_report-{args.job_id}.html",
        results / f"pipeline_timeline-{args.job_id}.html",
        results / f"pipeline_dag-{args.job_id}.html",
    ]
    html_error_tokens = ("traceback (most recent call last)", "internal server error", "cannot read properties of undefined")
    for html_path in current_html:
        if not html_path.is_file():
            check(f"Current execution HTML exists: {html_path.name}", False, str(html_path))
            continue
        html_text = html_path.read_text(encoding="utf-8", errors="replace").lower()
        detected = [token for token in html_error_tokens if token in html_text]
        check(f"Current execution HTML clean: {html_path.name}", not detected, f"error_signatures={detected}")

    # Run or reuse the exact-container exhaustive deep audit.
    deep_out = output_root / f"PGTK-deep-audit-{args.job_id}"
    if not args.reuse_deep_audit:
        audit_script = output / "embedded_audit_pgtk_results_parallel.py"
        audit_bytes = zlib.decompress(base64.b85decode(EMBEDDED_DEEP_AUDIT_B85.encode("ascii")))
        audit_digest = hashlib.sha256(audit_bytes).hexdigest()
        if audit_digest != EMBEDDED_DEEP_AUDIT_SHA256:
            check("Embedded deep-audit integrity", False, f"expected={EMBEDDED_DEEP_AUDIT_SHA256}, observed={audit_digest}")
        else:
            audit_script.write_bytes(audit_bytes)
            check("Embedded deep-audit integrity", True, audit_digest)
            audit_env = dict(os.environ)
            audit_env.update({"PGTK_SOURCE": str(project), "PGTK_RESULTS": str(results), "PGTK_AUDIT_OUT": str(deep_out), "PGTK_JOB_ID": args.job_id, "PGTK_MAX_EVENTS": str(args.max_events), "PGTK_AUDIT_WORKERS": str(workers), "PGTK_TARGET_IDS": "", "PGTK_STRICT": "false"})
            if args.pipeline_mode:
                deep = run([sys.executable, str(audit_script)], cwd=project, log=output / "deep_audit.log", env=audit_env)
            else:
                command = [str(args.apptainer), "exec", "--cleanenv", "--no-home"]
                for binding in dict.fromkeys([project, results, output, output_root]):
                    command.extend(["--bind", f"{binding}:{binding}"])
                command.extend([str(image), "env", *[f"{k}={v}" for k,v in audit_env.items() if k.startswith("PGTK_")], "python3", str(audit_script)])
                deep = run(command, cwd=project, log=output / "deep_audit.log")
            check("Exact-container exhaustive deep audit execution", deep.returncode == 0, f"exit={deep.returncode}")

    summary_path = deep_out / "summary.json"
    issues_path = deep_out / "issues.tsv"
    if summary_path.is_file() and issues_path.is_file():
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        deep_issues = read_tsv(issues_path)
        levels = Counter(row.get("level", "UNKNOWN") for row in deep_issues)
        metrics["deep_audit"] = summary
        check("Deep audit worker count", int(summary.get("audit_workers", 0)) == workers, f"reported={summary.get('audit_workers')}, requested={workers}")
        check("Deep audit errors", levels.get("ERROR", 0) == 0, json.dumps(dict(levels), sort_keys=True))
        check("All findings selected by deep audit", int(summary.get("display_events_selected", -1)) == int(summary.get("findings", -2)), f"selected={summary.get('display_events_selected')}, findings={summary.get('findings')}")
    else:
        check("Deep audit outputs", False, f"summary={summary_path.is_file()}, issues={issues_path.is_file()}")

    # Final report and archive.
    failures = [row for row in checks if row["Status"] == "FAIL"]
    warnings = [row for row in checks if row["Status"] == "WARN"]
    metrics.update({"elapsed_seconds": round(time.time() - started, 3), "checks": len(checks), "failures": len(failures), "warnings": len(warnings)})
    write_tsv(output / "checks.tsv", checks, ["Check", "Status", "Critical", "Details"])
    (output / "summary.json").write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    report = [f"# PGTK complete validation: job {args.job_id}", "", f"Overall: **{'PASS' if not failures else 'FAIL'}**", "", f"Workers requested and verified: **{workers}**", "", "## Checks", ""]
    report.extend(f"- **{row['Status']}** {row['Check']}: {row['Details']}" for row in checks)
    report += ["", "## Scope", "", f"- Declared Nextflow processes: {len(declared)}", f"- Observed process names: {len(observed)}", f"- Trace task rows: {len(trace_rows)}", f"- Samples: {', '.join(samples)}", f"- Validation elapsed seconds: {metrics['elapsed_seconds']}", "", "This validation proves source/runtime contracts and published-result integrity. It does not independently establish clinical validity or execute browser JavaScript."]
    (output / "REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    checksum = output / "checksums.sha256"
    with checksum.open("w", encoding="utf-8") as handle:
        for path in sorted(output.rglob("*")):
            if path.is_file() and path != checksum:
                handle.write(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(output)}\n")
    archive = output.with_suffix(".tar.gz")
    with tarfile.open(archive, "w:gz") as handle:
        handle.add(output, arcname=output.name)
    archive_sha = Path(str(archive) + ".sha256")
    archive_sha.write_text(f"{hashlib.sha256(archive.read_bytes()).hexdigest()}  {archive.name}\n", encoding="utf-8")

    print(f"OVERALL: {'PASS' if not failures else 'FAIL'}")
    print(f"FAILURES: {len(failures)}")
    print(f"WARNINGS: {len(warnings)}")
    print(f"REPORT: {output / 'REPORT.md'}")
    print(f"ARCHIVE: {archive}")
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
