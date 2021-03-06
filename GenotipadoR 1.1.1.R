

GenotipadoR=function(Ref="SARS-CoV-2.fasta", File="Query.fasta", Clades=read.table(file = 'clades.tsv', sep = '\t', header = TRUE), wsize, Crit=3, TestQuery, SiteQuery=TRUE)

## Funci�n con argumentos "Ref" (Nombre del archivo FASTA que se usar� como referencia), "File" (secuencia que se desea genotipar en formato FASTA), "Clades"
## (Archivo del que se leer�n los datos de genotipado) y wsize (Tama�o de ventana escogido). El par�metro "Crit" determinar� el criterio para
## considerar positivo un alineamiento. Se considerar� positivo aquel cuya puntuaci�n sea mayor a la Media + Crit * Desviaci�n Est�ndar.
## El argumento TestQuery permite introducir un vector de caracteres en lugar de un archivo .fasta.
## Los argumentos FullClades y SiteQuery permiten personalizar los resultados, eligiendo si ver o no todos los datos de los clados encontrados y si se desean ver las mutaciones respecto a posiciones gen�micas o de la secuencia introducida.

{

	require(seqinr)		## Comprobaci�n de los paquetes necesarios.

	require(Biostrings)

	Genoma=getSequence(read.fasta(Ref))[[1]]		## Define la variable "Genoma" como la secuencia extra�da del FASTA indicado.

	Query=getSequence(read.fasta(File))[[1]]

	if(!missing(TestQuery))

		{

		Query=TestQuery

		}

		if(missing(wsize))

		{

		wsize=length(Query)

		}

	AlineaQuery=c2s(Query)	## Almacena "Query" como secuencia de caracteres para el alineamiento.

	nventanas=(length(Genoma)%/% wsize) ## Divide el genoma en fragmentos de 100 nucle�tidos, sin incluir el resto de la divisi�n.

	resto=length(Genoma)%%wsize ## Almacena el resto de la divisi�n anterior.

	if (resto>0)

		{

		Ventanas <- vector(mode = "list", length = nventanas+1) 					## Crea una lista de longitud nventanas, y de
																## una posici�n m�s si la divisi�n no fue exacta,
		Ventanas[[nventanas+1]]=Genoma[((wsize*nventanas)+1):((wsize*nventanas)+resto)]	## en cuyo caso incluye en ella las posiciones
																## del resto de la divisi�n.
		 }							  				
									 				
	else

		{

		Ventanas <- vector(mode = "list", length = nventanas)

		}

	Ventanas[[1]]=Genoma[1:wsize]	## Crea la primera posici�n de la lista de fragmentos del genoma.

	for (i in 1:nventanas-1)      ## Bucle que genera el resto de las posiciones de la lista.

		{ 

		Ventanas[[i+1]]=Genoma[((wsize*i)+1):(wsize*(i+1))]

 		}

	DataFragmentos=data.frame	(						## Genera la data frame en la que se almacenar�n los resultados de los alineamientos y un color asignado seg�n su puntuaci�n.	

						"pAScore"=rep(0,length(Ventanas)),	## Columna de puntuaciones (vac��a de momento).

						"Color"=rep(NA,length(Ventanas))	## Columna de colores (vac��a de momento).

						)
	
	for (n in 1:length(Ventanas))	## Bucle para rellenar las puntuaciones.

		{

		DataFragmentos[n,1]=pairwiseAlignment(p=AlineaQuery,s=c2s(Ventanas[[n]]),type="global",scoreOnly=TRUE)

		}

	Mediana=median(DataFragmentos[,1])	## Almacena la media de las puntuaciones.

	RIQ=IQR(DataFragmentos[,1])		## Almacena el rango intercuart�lico de las puntuaciones.

	if (nrow(DataFragmentos)<10)		## Si el n�mero de ventanas es demasiado peque�o, modifica artificialmente la desviaci�n est�ndar para crear un criterio imposible de superar,
							## ya que igualmente en una muestra tan peque�a los resultados no ser�n muy informativos.
		{

			RIQ=Inf

		}

	Pos=Mediana+Crit*RIQ				## Busca positivos aplicando el criterio at�pico a la magnitud definida por el argumento "Crit". Por defecto se aplica el criterio at�pico leve ("Crit"=3).

	for (n in 1:length(Ventanas))			## Bucle que eval�a cada puntuaci�n y almacena un color seg�n si supera o no el criterio.

		{

		if (DataFragmentos[n,1]>Pos)

			{

			DataFragmentos[n,2]="Green"	## Si lo supera, almacena el color verde.

			}

		else

			{

			DataFragmentos[n,2]="Blue"	## Si no lo supera, almacena el color azul.

			}

		}
	
	if (all(DataFragmentos[,"Color"]!="Green"))	## Si ning�n fragmento supera el criterio, se considera positivo en su lugar el fragmento de puntuaci�n m�xima y emite un aviso.

		{

		DataFragmentos[which.max(DataFragmentos[,"pAScore"]),"Color"]="Yellow"
		
		message("Ninguna puntuaci�n del alineamiento supera el criterio establecido o el n�mero de ventanas es demasiado peque�o.
			Tomando como positivo en su lugar la puntuaci�n m�xima.
			Se recomienda ajustar el tama�o de ventana y/o el criterio.
				")

		}
	
	plot	(									## Construcci�n del gr�fico.

		1:length(Ventanas),DataFragmentos[,"pAScore"],

		main="Localizaci�n del fragmento",xlab="Ventana",ylab="Puntuaci�n",

		type="p", col=DataFragmentos[,"Color"]			

		)

	lines(DataFragmentos[,"pAScore"],col="black")			## Conecta los puntos con l��neas negras.

	CladeData=as.data.frame(Clades)[,-2] ## Carga como data frame los datos de los clados.

	VentanasPositivas=which(DataFragmentos[,"Color"]=="Green"|DataFragmentos[,"Color"]=="Yellow")	## Almacena la ventana cuyo alineamiento se ha considerado positivo.

	if(VentanasPositivas[1]==1&VentanasPositivas[length(VentanasPositivas)]!=length(Ventanas))	## Ampl�a con las ventanas adyacentes si el tama�o de ventana no coincidiera
																	## con la longitud de la secuencia introducida. Incluye ciertas condiciones
		{															## para no salirse del rango de ventanas.

		VentanasPositivas=c(VentanasPositivas,VentanasPositivas[length(VentanasPositivas)]+1)

		}

	if(VentanasPositivas[1]!=1&VentanasPositivas[length(VentanasPositivas)]==length(Ventanas))

		{
		
		VentanasPositivas=c(VentanasPositivas[1]-1,VentanasPositivas)

		}

	if(VentanasPositivas[1]!=1&VentanasPositivas[length(VentanasPositivas)]!=length(Ventanas))

		{
		
		VentanasPositivas=c(VentanasPositivas[1]-1,VentanasPositivas,VentanasPositivas+1)

		}

	SecuenciaPositiva=character(0)	## Se almacena la secuencia correspondiente a las ventanas seleccionadas.

	for (n in VentanasPositivas)

		{
	
			SecuenciaPositiva=c2s(c(SecuenciaPositiva,Ventanas[n][[1]]))

		}


	Alineamiento=pairwiseAlignment(p=AlineaQuery,s=SecuenciaPositiva,type="local",scoreOnly=FALSE)	## Alinea la secuencia introducida con la regi�n del genoma que ha salido positiva.
	
	QueryAlineado=toupper(s2c(as.character(alignedPattern(Alineamiento))))	## Almacenamiento del resultado del alineamiento como secuencias de caracteres.

	GenomaAlineado=toupper(s2c(as.character(alignedSubject(Alineamiento))))

	SumaResto=start(subject(Alineamiento))

	if(QueryAlineado[1]!=toupper(Query[1]))	## Ajuste por un defecto en la funci�n pairwiseAlignment() por el cual si la mutaci�n se encuentra en la primera posici�n esta se omite.

		{

		SumaResto=start(subject(Alineamiento))-1	## Actualiza manualmente la primera posici�n del genoma que alinea.

		QueryAlineado=c(toupper(Query[1]),QueryAlineado)	## Se añade el nucle�tido omitido al alineamiento, tanto en la secuencia query como en la referencia.

		GenomaAlineado=c(s2c(toupper(SecuenciaPositiva))[start(subject(Alineamiento))-1],GenomaAlineado)

		}

	if(is.na(QueryAlineado[length(Query)]))	## Ajuste por otro defecto equivalente si la mutaci�n est� en la �ltima posici�n.

		{

		QueryAlineado=c(QueryAlineado,toupper(Query[length(Query)]))

		GenomaAlineado=c(GenomaAlineado,s2c(toupper(SecuenciaPositiva))[end(subject(Alineamiento))+1])

		}

	if(all(QueryAlineado==GenomaAlineado))		## Si la secuencia introducida alinea a la perfecci�n con el genoma de referencia, se detiene la funci�n al no ser necesario el genotipado.

		{

			stop("La secuencia introducida es id�ntica al genoma de referencia.")

		}

	PosicionesDistintasAlineamiento=which(QueryAlineado!=GenomaAlineado)	## Almacena las posiciones del alineamiento en las que hay diferencias.

	PrimeraPosici�nGenoma=((VentanasPositivas[1]-1)*wsize+SumaResto)

	PosicionesDistintas=PosicionesDistintasAlineamiento+PrimeraPosici�nGenoma-1	## Peque�o ajuste necesario por contarse las mutaciones a partir del 1 y no del 0.

	Mutaciones=QueryAlineado[PosicionesDistintasAlineamiento]	## Almacena las mutaciones.

	Genotipado=data.frame(site=PosicionesDistintas,alt=Mutaciones)	## Se crea un data frame que almacena las posiciones y la mutaci�n en dicha posici�n.

	CladePos=integer(0)	## B�squeda en los datos de clados proporcionados. Se almacenan las filas cuyas mutaciones y posiciones coinciden con las del alineamiento.

	Found=integer(0) 							## Se define un vector de iguales caracter��sticas en el que se almacenar�n las posiciones de las mutaciones genotipadas.

	for (n in 1:nrow(Genotipado))		## Bucle donde se comparan las filas del genotipado con las de la base de datos.

		{

			donde = which(CladeData[,2]==Genotipado[n,1] & CladeData[,3]==Genotipado[n,2])	## Almacena las filas que coincidan.

  			if (length(donde>0)) { CladePos = c(CladePos,donde); Found = c(Found,n) }	## Si se encuentran coindicendias, se almacena en CladePos el n�mero de fila de la base de datos
																## donde se ha encontrado. Tambi�n se almacena un vector con las posiciones en la data frame
																## generada por la funci�n.
		}


	if (length(Found)>0)	## Si hay coincidencias entre el genotipado y la base de datos, se retiran de la data frame de genotipado.

		{

			Genotipado=Genotipado[-Found,]

		}

	if(nrow(Genotipado)>0)	## Si tras retirar las coindicencias a�n quedan elementos en la data frame de genotipado, se les asignan clados desconocidos.

		{

			Desc=data.frame(clade=rep("?",nrow(Genotipado)))

			Genotipado=cbind(Desc,Genotipado)

			message("Se han encontrado mutaciones no incluidas en la base de datos introducida. Se recomienda comprobar si est� actualizada.
					")

		}

		FoundClades=CladeData[CladePos,]

		Resultado=rbind(FoundClades,Genotipado)	## Almacena como resultado la informaci�n correspondiente a las mutaciones encontradas.

		if(SiteQuery)					## Si el argumento SiteQuery tiene el valor de TRUE, se actualizan las posiciones de las mutaciones para reflejar la posici�n en la
									## secuencia genotipada.
			{

			Resultado[,"site"]=(Resultado[,"site"]-PrimeraPosici�nGenoma)+1

			}


		message("Resultados del genotipado:
				")					## Muestra los resultados del genotipado.

		return(Resultado)
	
}