df = CSV.read(joinpath("data", "hcv.csv"), DataFrame)

stan_data = (;
    nt=length(df.time),
    time=df.time,
    cmt=coalesce.(df.cmt, 0) .+ 1,
    amt=coalesce.(df.amt, 0),
    rate=coalesce.(df.rate, 0),
    evid=df.evid,
    addl=zeros(Int64, nrow(df)),
    ss=df.ss,
    ii=coalesce.(df.ii, 0),
    nObs=nrow(filter(row -> row.evid == 0, df)),
    iObs=(@chain df begin
        @transform :iObs = 1:nrow(df)
        @rsubset :evid == 0
        @select :iObs
    end).iObs,
    nSubjects=length(unique(df.id)),
    step=15,
    start=[1, 1 + 15, 1 + 15 + 15],
    finish=[15, 15 + 15, 15 + 15 + 15],
    yPK=dropmissing(df, :yPK).yPK,
    yPD=dropmissing(df, :yPD).yPD
)

stan_init = (;
    logthetaKa=-0.2231435513142097,
    logthetaKe=-1.8971199848858813,
    logthetaVd=4.605170185988092,
    logthetan=0.6931471805599453,
    logthetaδ=-1.6094379124341003,
    logthetac=1.9459101490553132,
    logthetaEC50=-2.120263536200091,
    omegaKa=0.776261792561888,
    omegaKe=0.776261792561888,
    omegaVd=0.776261792561888,
    omegan=0.776261792561888,
    omegaδ=0.776261792561888,
    omegac=0.776261792561888,
    omegaEC50=0.776261792561888,
    sigmaPK=0.6895956605642456,
    sigmaPD=0.6895956605642456
)