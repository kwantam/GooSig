#!/usr/bin/python
#
# (C) 2018 Riad S. Wahby <rsw@cs.stanford.edu>

import hashlib

def _gen_group_obj(n, g, h):
    return type('',(object,),{"modulus": n, "g": g, "h": h})()

class Defs(object):
    winsize = 6     # for RSA operations
    combsize = 14   # for RSA operations
    hashfn = hashlib.sha256
    chalbits = 128  # V's challenge size

    # this is the list of primes from which P can choose to prove ability to take sqrt mod her RSA modulus
    # each one is a QR mod N with probability 1/2, so this list suffices except with vanishing probability
    primes = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109
             ,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233
             ,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367
             ,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499
             ,503,509,521,523,541,547,557,563,569,571,577,587,593,599,601,607,613,617,619,631,641,643
             ,647,653,659,661,673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797
             ,809,811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,941,947
             ,953,967,971,977,983,991,997]

    gen_group_obj = staticmethod(_gen_group_obj)
    # this modulus is the RSA-2048 challenge number
    Grsa2048 = _gen_group_obj(25195908475657893494027183240048398571429282126204032027777137836043662020707595556264018525880784406918290641249515082189298559149176184502808489120072844992687392807287776735971418347270261896375014971824691165077613379859095700097330459748808428401797429100642458691817195118746121515172654632282216869987549182422433637259085141865462043576798423387184774447920739934236584823824281198163815010674810451660377306056201619676256133844143603833904414952634432190114657544454178424020924616515723350778707749817125772467962926386356373289912154831438167899885040445364023527381951378636564391212010397122822120720357, 2, 3)
    # this modulus is the RSA-617 (also 2048-bit) challenge number)
    Grsa617 = _gen_group_obj(22701801293785014193580405120204586741061235962766583907094021879215171483119139894870133091111044901683400949483846818299518041763507948922590774925466088171879259465921026597046700449819899096862039460017743094473811056991294128542891880855362707407670722593737772666973440977361243336397308051763091506836310795312607239520365290032105848839507981452307299417185715796297454995023505316040919859193718023307414880446217922800831766040938656344571034778553457121080530736394535923932651866030515041060966437313323672831539323500067937107541955437362433248361242525945868802353916766181532375855504886901432221349733, 2, 3)
    # this modulus is from the now-defunct "America Online Root Certification Authority 2" certificate
    Gaol = _gen_group_obj(833287539555655830295437478623068715732893519470999001078575227726720833974299793327589892999446410695482689087817542178001025518000856645824542761860050201456988012423476797325416284135543115147497848693422016716844588890303978917861768153893481982002390666216248428266868192981018727094419599486404703564746373249295690184670845475555790369580007733124536111677885986700228007632811605375037348867111246742728259987286763650744497554590557328181022841922299670291568009554416925270309393290765959364224688474771305284758562184978264357696428875457896007139183128831260273238163285901042204479134355972621131905136584620604662058736669132342594748834225518170435481937761517048000372876190271778371991934048522521335509931794499153333036189250300411277900706695834937476478054354004604558199648562246970757390191511318015222081517956425652654710134304173256536041609721137081358431402394710892282816253447201619670765451884361947176326116403134771384774651137867814106943339025901079928519650268247936103502045262497449460438120483693210976009001793748449954988928005941680356289900003413781980827653292214516002118874375945551854673323131286257033477656922596627315766194604750898234659893295727170509283506976786630106108531865703, 2, 3)

    # some random primes for testing (saves time vs generating on the fly)
    primes_1024 = [ 160104197413736169533482759075458845682436887326335627265951716533753696521217542066928770121419596131312569877508547873037068057643306417272906101678445492650136649349831406458729482165924107830600707254298690550622589262441989507791406389691677212415870455495626407554897654029460469218783676336130440662717
                  , 52828911461112478441270628238860061511086340665919846590064620799570780148469340944817277481178659426571021663858644754300762202642019137175572336082456511752263174197160867137875968843838492616465556555003697225833924447576372394029110849757285642577712619586019143079206838681507774558984122215664121534097
                  , 30817232815047642690236588208685313400154297235718280962353461088444039444767896432041584053877144776813350642179762766786200442462726417721398882012241658852450801783557114398190455668187615608420407603137227573537518972006447022723327819681064903298527596961298903632728016406327572510976474173587383508441
                  , 56274201999185253089085184360489614464587875600149878584244396627532609037349612090400754567365855244937901931566175395708309079701360682759565279697244503974698875900562080663185314141329470374192565937999938483985541046278291016570306508184902225765396481128862669276234521449975232892841665242434137701527
                  , 50234396392114601906369517490271405706779453387174004387887062135354397816073227664596117010731087706650982551180551742319768285697770907971884622981786899412377283615833797740867528687626853264257923194614045996817479013499341140009125067651594245366245571102783530284519583414762906878214562523881269085501
                  , 166026642875106375051034241518019827798440232677503751158252527796596266221483628854205778090754395709550370401757546587718077868745599449976792599466486738662702419673857207822750830069904057032427140245821423345247375275338399989102158265070361664184890653627020491663664701444313643989318739917517437739597
                  , 120429362240799794734097056121781111889416344846617491446072479493520317897247381775672807532312222866461352013649767497341684666187080854322745664256348123768976812604407009553267600258149639792624030786554001031971876793707887247956012683117418500433216543529498241794347902218768939144492131816414949790259
                  , 111264220151166897989689719414310608509475706129310699955279411834403734180887710205473894077261789057437173161031291574227762408936952460764742602612178582974691641277043186165163571150528422777864796412585291171794132268692186881308777059073551575344943350511375411975719556347168232986530955605075295395337
                  , 99497545659105903710372261880208221140666317320992291847386569571982434788164598641016617307908621100911117217831675986177793230188626246852856478132652458129822999179806820189312391402489045092076084315988854658329929233688717839661764092300845275320364611782131314533387913425474034296637397377282061034179
                  , 35919557373817058082032299551948299441296544220659794348467889288910634536267059101853019831126857410815227627750661273866727497548261337508572797549774507296989199436139515090326615912551861402679696715157117084973624372549166565796845592843813491655279107597701276852654795628230159107195052415523909677043
                  , 77207464062561793441862919975165882232393310680694163164453026921299299176207348298180615922267593986045582691397012342285298072335394635283845297930356017429566214196325916661638620412222370592511751327358978532506292784251220432104858938892410296648581766602218067583077494942412209140853644721561241182617
                  , 20976502932218302879958781051144060332373134257432147161821578606633264347804545941924456565833047942802104632109179477068242991635914390983830918250054251032136774042849843585576701818168222059868334397817932524939950233894741832949244189818495652070981520002109196165488795079890260194760111215232654093213
                  , 64792968451745325196061747473710129984883596869894123714595765294707977204287966360455774750884867159416441219255556687182510974031342036015442773280448409632048117178171233343789717735428741934709554474378820724648016469938198950712248225804518723124281768472922202356681229452605304858007616298005349571523
                  , 146345619690110583716099072965451485680708140237410107816966158735120453142417798978390256999756514198972414950772966204670030109172227305876754018949768747012321003925986003404397948538364890140981122767689566656075928535963939736825266942979310599873470528792919336316267862397569308203938790815911491583817
                  , 10698280100175128261355784512002372257803990757307528706975691477763419472588738973707249180554046968151427562156098787430798191626464327694376410155646169916276988330640670988699299287735298651425125982065104959337595363675167163164429381660089648013622014505387208980462338992036927276635250678637991900059
                  , 144727526348696960777687188957623344586715633236231280926210579744117307276173206229818866025865948526274565132711325190982518978405490461840295029011450679038012329034465851473063664377365267450191200591372341286072154829271807914168564708600945013830189972450612656059143580154493679883768964063169455150451
                  ]

    primes_2048 = [ 25847034201546780436474342804132278483696270887601946305801448178962948965302643193578313051153708440165586169446573111899545638276565225547370613926376202239986186394490817640433772010604479374729636877928340132492900126597741285585650417070732892578335460551215808711573626798728360522677447511895214001659196380148436581780985161666118027818353700314750304649102693217892523951160234688961921833921702015271727119994786036544348110728236000356810227691086994146439877828093383576478999943334939473373693932031913683006851237810106317060060670987149099666470644760334604393839712127685796326237126412086058517030651
                  , 16685694124116278119373951214115941534612189784654380837587853678710895800990248828598342503403539634780917184477426245693065888489763794871485405488028663760341682804666885555994273589405500760381869075719153999265208572379444102825725841709132945102716055341328506203013850453402399053957744562273061385819793837259385753736376678606576687959931908923712451329564505988408329613421459701883273638713598663301570848574671765575928619812340182346721265286950878078799706570215472403257353425479962496988595498280185125818449933136603190472852628969582278763374624868131590297556310302226188860754096190467807090040877
                  , 14435992843798391187135161818224243192437496549297003937503687162256087387658357332239013638177828388393296280174030574325182448030091558695212026230455574487126069320558229870034061565619621976653307666888412638034825630437564225397242431656561321584306654147689751275220761780421307459533854130263824229281621813714768129412655080414004541312303915634447187145820931279057684094503308571786051887579117705962966360644935130466960460008604005384237648774068571118696185055316634955685332389901129326246258081088199625692290271359195269026755421365230353530576952119906799097959359461177579370547311513459151723280309
                  , 18404707629122352924109127109314440313634545687157214353243891214755044740041921405008260565120219441482945462173737403425925642037560530985000533659799801216836842601441135207520665496053735858465971575613776158696648708229189937217681733567121747948035703833819725139693110676735881528922108150520638721743371521744752393830109996305682602665400043766045412805121692573326731522634021825305086771004483680194330317928227248177525822570109293880446950596763207130955071177582866116042337555976099654007985233308539137815986106042851717459644495744128668241306225455499441210832901024346912674577462918426353433486599
                  , 5280103075115904558508765027567682103999799385075455729255313607740458670654451837488656661174870040131125195302568330903150028901612311610861488437468541273699109061240613794257202167700224092506149782472810790069581159146555963757196008297274492405757021568945618872731291099830758072469400915268030365390441456640716106771499041296888399540384763822586997581833623839432859912972998901579800671303641738960486256160616767548343713500120554476891956567755013069645048733681623498719449738180003654383395349213770892680591175041711754900441630190345732428323095998987986352061632297615433328449749193392849707828769
                  , 22667104471494592288099945177752692694153971164710449416170916790713538459827167064217679971049921583125499186663152505913982536488228996181467965051766512759089976041483192994651459059238971155334217890518602198953607609141984523884974100067770827231818319378963775510379785999908902993219206649314791136060662744080989551770784503823001731700025141492153124022857989221728021434250112197366254827028502753749880291215209747013295418076184447599255402976615333143567703017984736301780852720468357267573591524310702340568397901845986378996145081540148880393248385356826379946298212061289362695190770872627441860286223
                  , 31369525209554606243781872919147324495102438235204029024209863960990364788173016485609984157542616821530687077860054182657226524268508850870308695525705464179120260058399486119923345168413129341961927332012856973037798959880697235803480704035779279173240855808158142157971123503154296708811317754274213461651382233302344648050509226964917724185500045556091237800513090875435392081653086625727098741948718259026787073395040413243577779586443920174270515835168155527724979839452996835774962902591815840849991142745534394197953906778975966142467432235494552657288854307685522009043975433801213110572831125772428981005613
                  , 9830889462965132138705959462261553388126280626805187143388562744640198373566296292273357048632589230670080127192008901370831775879317326960347961482833817171115437262853880339669403198417241594166879587414133075107451372037847631227713408882880398509149327046749725833772140201359658190669764465796820928318881602091170613212182342256385965652082031372762599856552019784468687370788837279870703159968394203891693401628970626854333045166966465593433831533032912244431989114902745938091527792128853799787836770246456315935650842345253375010353964349164641991917776480058108762281089941261153346857167120797701797173217
                  , 16425001985705049801406567854047086073127739400886126153156457090314190058792977582132326854862479031884437919286067612592750885780324499346040778117136571506414806055640673328850574898234716103123831451036219427744982197364896110018082074373560968723240016617417761480092490817533728673675894824760064888029030969786857427299891767724464100751425180717972298635050804633282834206479831362178462016634908982950618850270683421048435990333914799277540793592187097212465985935361520816006226206142522253660049411914334843188999273366975679945638123242981045376686086415351857926240087756000529444249797948281140860750797
                  , 12174000277072663132694111281367365428645877576618290061201151401411879625094186161225139001483127756850347306901565160901724456236850818840308128553980292264151202304387603539937849158327594379952903498414842500096909122346854920785464176200030890325516276678069322996265407928741941020271815742604393440830668186520437195654768441331465412571264473078339881537597636077932794077180222506207121343238963393278270853526704883650546969581986807346741158197537368810933762588613969529258234138213444363793559806753896940263386806544083417836806265688209717892609305307822002489081816938971540715010204635638959872069383
                  , 3237103031500867098321186786517528168574799789915832107232898654735840604943570992286132082345937243291531205211192012779075563641326889836973070545211822990617652558426452697857789899999279683993090012117464672117886770602440168319863143917822916206683351635354555063916015816370917232268413925159859287908508182052893407531427829957936920793714202088216692972386200887736250286344513457226319963902073528399561281693746379915787996867612570396498063107400269305088714514781363835891303118394643451756680125341238941200325994765342167708278057001805841012849015739791345634583206018308274378034775423507486682291441
                  , 23798395843759865970572157405719299409679832080031916162437000667138444865321976549957060252027933852435488114051971679262317660941838502637952049563627504489830541114255355055173575438534085749093838573154617770199255452691401501447448985499390415175837609475399831264415203639498953614365358450310849286047470753782408733509604190681375662715761198453929075911524651421748647975291134925207394558178441805235091431738981785099338589230380955599231105404180359080673242469902417933923785998557390169755852795155240452383094318092140045538401186275331485291236035719677943588391038583107742219935264275607723334588021
                  , 4180671455942636543094139596341252547581638085999390538168781596666986312916316791586350549270383786652688738133561023106628512993555465905685635185509792172540769848856186647779433181360435282821345145176892498810050194148836627533777703015698219288642887929879987091108360342564013772571067705932172721125976953679624100542018155417064037394044882992756154197852850976377939745237288587249405604225705042846401919591830021713507244051120288238344392109854312502233677349574017528999827121159888999209698851606749628795460677993125630927672818963553720538572517969079138105020363675169913748497505438535096129094543
                  , 13498311400665465113692361968715863060943526907836574816844702975139163986008331077392886710016607552225961987734904867853270635328516523541804454220536047246050744924864205915177960153489189762088771785771792330724995531512702565171713426397860490262184371822873702155319146672894031125790600394666890016397666349827660063420547801489067037013417986070508442887833624986374625894727010690822525602414452704243663855059608392832571057586003290565306256647718766131862247517539906186862423783643549573629253881416830342642043364533411602933262834266547002676886644504668871898870606289473798362201631868550391485261817
                  , 23327344558654898417807202127235596756043853075449217195800374285863570109994985422746315115613399388482819306231463673042914944001635008012441817938140619119975695355879792660806335386219100733604935187312996604320514233616080226114662990106230206895848550953635979359397598700137432346931335375082648318008317432443769746291452408552230469862815593503779549612181045113917334290858886861384958798931545668122178627783564003535532440937362056754015627420708871620538761982919958787428669649550623931592001608149255309112723581989853951586189622828078638304544957659516956904162247400293872782370874168557361801709499
                  , 7797475962205081870379191226238756363914735330638510142979001385759278289428172500433219492402719413813213353794634519914085596871088029020930597102155050063349311475963159312523950132864740760938472382180770869970444905454393697268441434535356975507854658375869242166191654801621568857412166088953956766961530000830007641956183444570307155650450051979464775690617245426905450694525265441731687586290221525296981187118531950369082343812021032829166615628759440051417737697463983449619956441795607846427445827011773838819989877116539203848895055939800852004184549768172936393304644932490210448159285118834217596434563
                  ]