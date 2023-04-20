#include "HSS.h"

void Paillier_Gen(PK &pk, ZZ &d)
{
    ZZ p, q, Phi_N, temp1, temp2, h;
    
    // GenGermainPrime(p, 1536); // safe prime
    p = conv<ZZ>("2410312426921032588580116606028314112912093247945688951359675039065257391591803200669085024107346049663448766280888004787862416978794958324969612987890774651455213339381625224770782077917681499676845543137387820057597345857904599109461387122099507964997815641342300677629473355281617428411794163967785870370368969109221591943054232011562758450080579587850900993714892283476646631181515063804873375182260506246992837898705971012525843324401232986857004760339321639");
    // GenGermainPrime(q, 1536);
    q = conv<ZZ>("2410312426921032588580116606028314112912093247945688951359675039065257391591803200669085024107346049663448766280888004787862416978794958324969612987890774651455213339381625224770782077917681499676845543137387820057597345857904599109461387122099507964997815641342300677629473355281617428411794163967785870370368969109221591943054232011562758450080579587850900993714892283476646631181515063804873375182260506246992837898705971012525843324401232986857004760339319223");
    
    mul(pk.N, p, q);        // N = p * q
    mul(pk.N2, pk.N, pk.N); // N^2

    sub(temp1, p, 1);         // phi(p)
    sub(temp2, q, 1);         // phi(q)
    mul(Phi_N, temp1, temp2); // phi(N)

    InvMod(temp1, Phi_N, pk.N);
    mul(temp1, temp1, Phi_N);
    mul(temp2, pk.N, Phi_N);
    rem(d, temp1, temp2);

    pk.DJN_OPEN = false;
    pk.Pre_OPEN = false;
    if (pk.DJN_OPEN)
    {
        RandomBnd(h, pk.N);
        MulMod(h, -h, h, pk.N2);
        PowerMod(pk.hs, h, pk.N, pk.N2);
    }
}

void Pailler_Enc(ZZ &ct, PK pk, ZZ x)
{
    ZZ r, temp1, temp2;
    if (pk.DJN_OPEN)
    {
        if (pk.Pre_OPEN) {
            temp1 = conv<ZZ>("90505564256596579447266009133768049721586049409607196041549895898546110759002030713205553200083864955831961581866107550059857274096253961912893165433885014381529808569987280315251114614015865868774181741051566703782611342175378062265789461624884596255507610917276421472800808710453396107856356503694272381040446824216432529278783262944100357547538971341785194991597273086477349386490143917418688972838451809351360096191914020620277989839746377345945341776373672638453452785810946598572062480058169059782337235752060260947225360714948973950140102922377233007015642616373667030659415840004260981588199911273839815793590359352885891837937380415810129282963125866776106874500314815666537088024436135690771542402638449005948616529757594621511079139137043850275052109466540602716120674851635540557239611050112172558827178076476962815982816494514778727941719897155951598168103647810743708861253495043981411055179945767385164387079800813922727476276167496744665746739319800029206034467538276907274184250739384556096490217750957929899995388618939724078479555156090671865648104093779426500688154764389258177790596416349250815436225815352897154896287026378550933268442553216007949657071389560602781325474209540973879658923771330063739110814562490052041279265269259993884832892123353797265917435756350386619388993457144029561917040593027198431215331955822939724507074592446011551005945741493136481769239059614903250953790195054396608459360504634466899120837837017295361439029382138854507384271482014654859089091344870570496589571837654475366008538191382503655001835719927537241978224693925186120131493323304842671420779221149419027399224268590751650882474154959683292966038795434525442588731101284563738267676119304361397712192278973066027449908483544309652697453005132293907031254185833418062380594466593333562623402850659250895796602931247918254228180579281229050556425659657944726600913376804972158604940960719604154989589854611075900203071320555320008386495583196158186610755005985727409625396191289316543388501438152980856998728031525111461401586586877418174105156670378261134217537806226578946162488459625550761091727642147280080871045339610785635650369427238104044682421643252927878326294410035754753897134178519499159727308647734938649014391741868897283845180935136009619191402062027798983974637734594534177637367263845345278581094659857206248005816905978233723575206026094722536071494897395014010292237723300701564261637366703065941584000426098158819991127383981579359035935288589183793738041581012928296312586677610687450031481566653708802443613569077154240263844900594861652975759462151107913913704385027505210946654060271612067485163554055723961105011217255882717807647696281598281649451477872794171989715595159816810364781074370886125349504398141105517994576738516438707980081392272747627616749674466574673931980002920603446753827690727418425073938455609649021775095792989999538861893972407847955515609067186564810409377942650068815476438925817779059641634925081543622581535289715489628702637855093326844255321600794965707138956060278132547420954097387965892377133006373911081456249005204127926526925999388483289212335379726591743575635038661938899345714402956191704059302719843121533195582293972450707459244601155100594574149313648176923905961490325095379019505439660845936050463446689912083783701729536143902938213885450738427148201465485908909134487057049658957183765447536600853819138250365500183571992753724197822469392518612013149332330484267142077922114941902739922426859075165088247415495968329296603879543452544258873110128456373826767611930436139771219227897306602744990848354430965269745300513229390703125418583341806238059446659333356262340285065925089579660293124791825422818057928122");
        }
        else {
            ZZ r_bnd = conv<ZZ>("2410312426921032588580116606028314112912093247945688951359675039065257391591803200669085024107346049663448766280888004787862416978794958324969612987890774651455213339381625224770782077917681499676845543137387820057597345857904599109461387122099507964997815641342300677629473355281617428411794163967785870370368969109221591943054232011562758450080579587850900993714892283476646631181515063804873375182260506246992837898705971012525843324401232986857004760339316736");
            RandomBnd(r, r_bnd);
            PowerMod(temp1, pk.hs, r_bnd, pk.N2);
        }
    }
    else
    {
        RandomBnd(r, pk.N2);
        PowerMod(temp1, r, pk.N, pk.N2);
    }
    //PowerMod(temp2, pk.N+1,x,pk.N2);
    MulMod(temp2, pk.N, x, pk.N2);
    add(temp2, temp2, 1);
    MulMod(ct, temp1, temp2, pk.N2);

}

void Pailler_Dec(ZZ &x, PK pk, ZZ sk, ZZ ct)
{
    ZZ temp;
    PowerMod(temp, ct, sk, pk.N2);
    sub(temp, temp, 1);
    div(x, temp, pk.N);
}

void HSS_Gen(PK &pk, EK &ek0, EK &ek1)
{
    ZZ d, temp1, temp2;
    Vec<ZZ> d_Bsk;

    Paillier_Gen(pk, d);

    pk.k = 1024;
    power(pk.Bmsg, 2, pk.k);
    // power(pk.Bsk, 2, pk.k);
    pk.Bsk = pk.Bmsg;

    temp1 = d;
    while (temp1 != 0)
    {
        // temp2 = temp1 % Bsk; temp1 = temp1 / Bsk
        DivRem(temp1, temp2, temp1, pk.Bsk);
        d_Bsk.append(temp2);
    }
    pk.l = d_Bsk.length();

    mul(temp1, pk.Bmsg, pk.Bsk);
    pk.D.SetLength(pk.l);
    int i;
    for (i = 0; i < pk.l; ++i)
    {
        RandomBnd(temp2, temp1);
        ek1.d_Bsk.append(temp2);
        sub(temp2, temp2, d_Bsk[i]);
        ek0.d_Bsk.append(temp2);
        Pailler_Enc(pk.D[i], pk, d_Bsk[i]);
    }
}

void HSS_Free(PK &pk, EK &ek0, EK &ek1)
{
    pk.D.kill();
    ek0.d_Bsk.kill();
    ek1.d_Bsk.kill();
}

void HSS_Input(Vec<ZZ> &Ix, PK pk, ZZ x)
{
    if (!Ix.length())
    {
        Ix.SetLength(pk.l + 1);
    }
    ZZ r, temp1, temp2;
    Pailler_Enc(Ix[0], pk, x);

    int i;
    for (i = 0; i < pk.l; ++i)
    {
        if (pk.DJN_OPEN)
        {
            if (pk.Pre_OPEN) {
                temp1 = conv<ZZ>("90505564256596579447266009133768049721586049409607196041549895898546110759002030713205553200083864955831961581866107550059857274096253961912893165433885014381529808569987280315251114614015865868774181741051566703782611342175378062265789461624884596255507610917276421472800808710453396107856356503694272381040446824216432529278783262944100357547538971341785194991597273086477349386490143917418688972838451809351360096191914020620277989839746377345945341776373672638453452785810946598572062480058169059782337235752060260947225360714948973950140102922377233007015642616373667030659415840004260981588199911273839815793590359352885891837937380415810129282963125866776106874500314815666537088024436135690771542402638449005948616529757594621511079139137043850275052109466540602716120674851635540557239611050112172558827178076476962815982816494514778727941719897155951598168103647810743708861253495043981411055179945767385164387079800813922727476276167496744665746739319800029206034467538276907274184250739384556096490217750957929899995388618939724078479555156090671865648104093779426500688154764389258177790596416349250815436225815352897154896287026378550933268442553216007949657071389560602781325474209540973879658923771330063739110814562490052041279265269259993884832892123353797265917435756350386619388993457144029561917040593027198431215331955822939724507074592446011551005945741493136481769239059614903250953790195054396608459360504634466899120837837017295361439029382138854507384271482014654859089091344870570496589571837654475366008538191382503655001835719927537241978224693925186120131493323304842671420779221149419027399224268590751650882474154959683292966038795434525442588731101284563738267676119304361397712192278973066027449908483544309652697453005132293907031254185833418062380594466593333562623402850659250895796602931247918254228180579281229050556425659657944726600913376804972158604940960719604154989589854611075900203071320555320008386495583196158186610755005985727409625396191289316543388501438152980856998728031525111461401586586877418174105156670378261134217537806226578946162488459625550761091727642147280080871045339610785635650369427238104044682421643252927878326294410035754753897134178519499159727308647734938649014391741868897283845180935136009619191402062027798983974637734594534177637367263845345278581094659857206248005816905978233723575206026094722536071494897395014010292237723300701564261637366703065941584000426098158819991127383981579359035935288589183793738041581012928296312586677610687450031481566653708802443613569077154240263844900594861652975759462151107913913704385027505210946654060271612067485163554055723961105011217255882717807647696281598281649451477872794171989715595159816810364781074370886125349504398141105517994576738516438707980081392272747627616749674466574673931980002920603446753827690727418425073938455609649021775095792989999538861893972407847955515609067186564810409377942650068815476438925817779059641634925081543622581535289715489628702637855093326844255321600794965707138956060278132547420954097387965892377133006373911081456249005204127926526925999388483289212335379726591743575635038661938899345714402956191704059302719843121533195582293972450707459244601155100594574149313648176923905961490325095379019505439660845936050463446689912083783701729536143902938213885450738427148201465485908909134487057049658957183765447536600853819138250365500183571992753724197822469392518612013149332330484267142077922114941902739922426859075165088247415495968329296603879543452544258873110128456373826767611930436139771219227897306602744990848354430965269745300513229390703125418583341806238059446659333356262340285065925089579660293124791825422818057928122");
            }
            else {
                ZZ r_bnd = conv<ZZ>("2410312426921032588580116606028314112912093247945688951359675039065257391591803200669085024107346049663448766280888004787862416978794958324969612987890774651455213339381625224770782077917681499676845543137387820057597345857904599109461387122099507964997815641342300677629473355281617428411794163967785870370368969109221591943054232011562758450080579587850900993714892283476646631181515063804873375182260506246992837898705971012525843324401232986857004760339316736");
                RandomBnd(r, r_bnd);
                PowerMod(temp1, pk.hs, r_bnd, pk.N2);
            }
        }
        else
        {
            RandomBnd(r, pk.N2);
            PowerMod(temp1, r, pk.N, pk.N2);
        }
        PowerMod(temp2, pk.D[i], x, pk.N2);
        MulMod(Ix[i + 1], temp1, temp2, pk.N2);
    }
}

void HSS_ConvertInput(Vec<ZZ> &Mx, int b, PK pk, EK ek, Vec<ZZ> Ix, int &prf_key)
{
    ZZ temp1;
    temp1 = PRF_ZZ(prf_key++, pk.Bmsg);
    Vec<ZZ> M1;

    if (b == 0)
    {
        M1.append(temp1);
    }
    else
    {
        AddMod(temp1, temp1, ZZ(1), pk.N);
        M1.append(temp1);
    }
    M1.append(ek.d_Bsk);
    HSS_Mul(Mx, pk, ek, Ix, M1, prf_key);
}

void HSS_Mul(Vec<ZZ> &Mz, PK pk, EK ek, Vec<ZZ> Ix, Vec<ZZ> My, int &prf_key)
{
    if (!Mz.length())
    {
        Mz.SetLength(pk.l + 1);
    }

    ZZ temp1, temp2, temp3;
    temp2 = 0;
    int i;
    for (i = 0; i < pk.l; ++i)
    { 
        power(temp1, pk.Bsk, i);  
        mul(temp1, temp1, My[i + 1]);
        add(temp2, temp2, temp1);
    }

    for (i = 0; i < pk.l + 1; ++i)
    {
        PowerMod(temp1, Ix[i], temp2, pk.N2);
        HSS_DDLog(Mz[i], pk, temp1);
        AddMod(Mz[i], Mz[i], PRF_ZZ(prf_key++, pk.N), pk.N);
    }
}

void HSS_DDLog(ZZ &z, PK pk, ZZ g)
{
    ZZ h1, h, temp1;
    DivRem(h1, h, g, pk.N); // h = g % N; h1 = g / N
    InvMod(temp1, h, pk.N);
    MulMod(z, h1, temp1, pk.N);
}

void HSS_AddMemory(Vec<ZZ> &Mz, PK pk, Vec<ZZ> Mx, Vec<ZZ> My)
{
    if (!Mz.length())
    {
        Mz.SetLength(pk.l + 1);
    }
    int i;
    for (i = 0; i < pk.l + 1; ++i)
    {
        add(Mz[i], Mx[i], My[i]);
    }
}

void HSS_AddInput(Vec<ZZ> &Iz, PK pk, Vec<ZZ> Ix, Vec<ZZ> Iy)
{
    if (!Iz.length())
    {
        Iz.SetLength(pk.l + 1);
    }

    int i;
    for (i = 0; i < pk.l + 1; ++i)
    {
        MulMod(Iz[i], Ix[i], Iy[i], pk.N2);
    }
}

void HSS_Evaluate(Vec<ZZ> &y_b_res, int b, Mat<ZZ> Ix, PK pk, EK ekb, int &prf_key, vector<vector<int>> F_TEST)
{
    y_b_res.kill();
    y_b_res.SetLength(pk.l + 1);
    Vec<ZZ> M1;
    ZZ temp1;
    temp1 = PRF_ZZ(prf_key++, pk.Bmsg);

    if (b == 0)
    {
        M1.append(temp1);
    }
    else
    {
        add(temp1,temp1,1);
        M1.append(temp1);
    }
    M1.append(ekb.d_Bsk);
    Vec<ZZ> Monomial;
    int i, j, k;
    for (i = 0; i < F_TEST.size(); ++i)
    {
        Monomial = M1;
        for (j = 0; j < Ix.NumRows(); ++j)
        {
            for (k = 0; k < F_TEST[i][j]; ++k)
            {
                HSS_Mul(Monomial, pk, ekb, Ix[j], Monomial, prf_key);
            }
        }
        add(y_b_res, y_b_res, Monomial);
    }
}