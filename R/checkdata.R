checkdata<-function(epg,database,typed.gts=NULL,reference.profiles=NULL)
{
fail<-FALSE
mdb<-unique(database[['marker']])
mepg<-unique(epg[['marker']])
if(!all(mepg%in%mdb)) 
	{
	warning('                    marker(s) used in epg but missing in database: ',
		paste(setdiff(mepg,mdb),collapse=' '))
	fail<-TRUE
	}
for(m in mepg) if(m%in%mdb)
	{
	miss<-setdiff(subset(epg,epg[['marker']]==m)$allele,subset(database,database[['marker']]==m)$allele)
	if(0!=length(miss)) 
	{
	warning('                    allele(s) used in epg but missing in database: ',miss,
		' on marker ',m,' ')
	fail<-TRUE
	}
	}
if(!is.null(reference.profiles))
	{
	mrp<-unique(reference.profiles[['marker']])
	if(!all(mepg%in%mrp))
		{
		warning('          marker(s) used in epg but missing in reference.profiles: ',
			paste(setdiff(mdb,mrp),collapse=' '))
		fail<-TRUE
		}
	for(m in mrp) if(m%in%mdb)
		{
		miss<-setdiff(subset(reference.profiles,reference.profiles[['marker']]==m)$allele,subset(database,database[['marker']]==m)$allele)
		if(0!=length(miss)) 
		{
		warning('     allele(s) used in reference.profiles but missing in database: ',miss,' on marker ',m)
		fail<-TRUE
		}
		}
	}
if(!is.null(typed.gts)) for(i in seq_along(typed.gts))
	{
	mtp<-typed.gts[[i]][['marker']]
	if(!all(mepg%in%mtp)) 
		{
		warning('              marker(s) used in epg but missing in typed.gts[[',i,']]: ',
			paste(setdiff(mepg,mtp),collapse=' '))
		fail<-TRUE
		}
	for(m in mtp) if(m%in%mdb)
		{
		miss<-setdiff(subset(typed.gts[[i]],typed.gts[[i]][['marker']]==m)$allele1,subset(database,database[['marker']]==m)$allele)
		miss<-union(miss,
			setdiff(subset(typed.gts[[i]],typed.gts[[i]][['marker']]==m)$allele2,subset(database,database[['marker']]==m)$allele))
		if(0!=length(miss)) 
		{
		warning('     allele(s) used in typed.gts[[',i,']] but missing in database: ',
			paste(miss,collapse=' '),' on marker ',m)
		fail<-TRUE
		}
		}
	}
if(fail) stop('missing markers or alleles')
}


