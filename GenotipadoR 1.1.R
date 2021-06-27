

GenotipadoR=function(Ref="SARS-CoV-2.fasta", File="Query.fasta", Clades=read.table(file = 'clades.tsv', sep = '\t', header = TRUE), wsize, Crit=3, TestQuery, FullClades=FALSE, SiteQuery=TRUE)

## Función con argumentos "Ref" (Nombre del archivo FASTA que se usará como referencia), "File" (secuencia que se desea genotipar en formato FASTA), "Clades"
## (Archivo del que se leerán los datos de genotipado) y wsize (Tamaño de ventana escogido). El parámetro "Crit" determinará el criterio para
## considerar positivo un alineamiento. Se considerará positivo aquel cuya puntuación sea mayor a la Media + Crit * Desviación Estándar.
## El argumento TestQuery permite introducir un vector de caracteres en lugar de un archivo .fasta.
## Los argumentos FullClades y SiteQuery permiten personalizar los resultados, eligiendo si ver o no todos los datos de los clados encontrados y si se desean ver las mutaciones respecto a posiciones genómicas o de la secuencia introducida.

{

	require(seqinr)		## Comprobación de los paquetes necesarios.

	require(Biostrings)

	Genoma=getSequence(read.fasta(Ref))[[1]]		## Define la variable "Genoma" como la secuencia extraída del FASTA indicado.

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

	nventanas=(length(Genoma)%/% wsize) ## Divide el genoma en fragmentos de 100 nucleótidos, sin incluir el resto de la división.

	resto=length(Genoma)%%wsize ## Almacena el resto de la división anterior.

	if (resto>0)

		{

		Ventanas <- vector(mode = "list", length = nventanas+1) 					## Crea una lista de longitud nventanas, y de
																## una posición más si la división no fue exacta,
		Ventanas[[nventanas+1]]=Genoma[((wsize*nventanas)+1):((wsize*nventanas)+resto)]	## en cuyo caso incluye en ella las posiciones
																## del resto de la división.
		 }							  				
									 				
	else

		{

		Ventanas <- vector(mode = "list", length = nventanas)

		}

	Ventanas[[1]]=Genoma[1:wsize]	## Crea la primera posición de la lista de fragmentos del genoma.

	for (i in 1:nventanas-1)      ## Bucle que genera el resto de las posiciones de la lista.

		{ 

		Ventanas[[i+1]]=Genoma[((wsize*i)+1):(wsize*(i+1))]

 		}

	DataFragmentos=data.frame	(						## Genera la data frame en la que se almacenarán los resultados de los alineamientos y un color asignado según su puntuación.	

						"pAScore"=rep(0,length(Ventanas)),	## Columna de puntuaciones (vací­a de momento).

						"Color"=rep(NA,length(Ventanas))	## Columna de colores (vací­a de momento).

						)
	
	for (n in 1:length(Ventanas))	## Bucle para rellenar las puntuaciones.

		{

		DataFragmentos[n,1]=pairwiseAlignment(p=AlineaQuery,s=c2s(Ventanas[[n]]),type="global",scoreOnly=TRUE)

		}

	Mediana=median(DataFragmentos[,1])	## Almacena la media de las puntuaciones.

	RIQ=IQR(DataFragmentos[,1])		## Almacena el rango intercuartílico de las puntuaciones.

	if (nrow(DataFragmentos)<10)		## Si el número de ventanas es demasiado pequeño, modifica artificialmente la desviación estándar para crear un criterio imposible de superar,
							## ya que igualmente en una muestra tan pequeña los resultados no serán muy informativos.
		{

			RIQ=Inf

		}

	Pos=Mediana+Crit*RIQ				## Busca positivos aplicando el criterio atípico a la magnitud definida por el argumento "Crit". Por defecto se aplica el criterio atípico leve ("Crit"=3).

	for (n in 1:length(Ventanas))			## Bucle que evalúa cada puntuación y almacena un color según si supera o no el criterio.

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
	
	if (all(DataFragmentos[,"Color"]!="Green"))	## Si ningún fragmento supera el criterio, se considera positivo en su lugar el fragmento de puntuación máxima y emite un aviso.

		{

		DataFragmentos[which.max(DataFragmentos[,"pAScore"]),"Color"]="Yellow"
		
		message("Ninguna puntuación del alineamiento supera el criterio establecido o el número de ventanas es demasiado pequeño.
			Tomando como positivo en su lugar la puntuación máxima.
			Se recomienda ajustar el tamaño de ventana y/o el criterio.
				")

		}
	
	plot	(									## Construcción del gráfico.

		1:length(Ventanas),DataFragmentos[,"pAScore"],

		main="Localización del fragmento",xlab="Ventana",ylab="Puntuación",

		type="p", col=DataFragmentos[,"Color"]			

		)

	lines(DataFragmentos[,"pAScore"],col="black")			## Conecta los puntos con lí­neas negras.

	CladeData=as.data.frame(Clades)[,-2] ## Carga como data frame los datos de los clados.

	VentanasPositivas=which(DataFragmentos[,"Color"]=="Green"|DataFragmentos[,"Color"]=="Yellow")	## Almacena la ventana cuyo alineamiento se ha considerado positivo.

	if(VentanasPositivas[1]==1&VentanasPositivas[length(VentanasPositivas)]!=length(Ventanas))	## Amplía con las ventanas adyacentes si el tamaño de ventana no coincidiera
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


	Alineamiento=pairwiseAlignment(p=AlineaQuery,s=SecuenciaPositiva,type="local",scoreOnly=FALSE)	## Alinea la secuencia introducida con la región del genoma que ha salido positiva.
	
	QueryAlineado=toupper(s2c(as.character(alignedPattern(Alineamiento))))	## Almacenamiento del resultado del alineamiento como secuencias de caracteres.

	GenomaAlineado=toupper(s2c(as.character(alignedSubject(Alineamiento))))

	SumaResto=start(subject(Alineamiento))

	if(QueryAlineado[1]!=toupper(Query[1]))	## Ajuste por un defecto en la función pairwiseAlignment() por el cual si la mutación se encuentra en la primera posición esta se omite.

		{

		SumaResto=start(subject(Alineamiento))-1	## Actualiza manualmente la primera posición del genoma que alinea.

		QueryAlineado=c(toupper(Query[1]),QueryAlineado)	## Se aÃ±ade el nucleótido omitido al alineamiento, tanto en la secuencia query como en la referencia.

		GenomaAlineado=c(s2c(toupper(SecuenciaPositiva))[start(subject(Alineamiento))-1],GenomaAlineado)

		}

	if(is.na(QueryAlineado[length(Query)]))	## Ajuste por otro defecto equivalente si la mutación está en la última posición.

		{

		QueryAlineado=c(QueryAlineado,toupper(Query[length(Query)]))

		GenomaAlineado=c(GenomaAlineado,s2c(toupper(SecuenciaPositiva))[end(subject(Alineamiento))+1])

		}

	if(all(QueryAlineado==GenomaAlineado))		## Si la secuencia introducida alinea a la perfección con el genoma de referencia, se detiene la función al no ser necesario el genotipado.

		{

			stop("La secuencia introducida es idéntica al genoma de referencia.")

		}

	PosicionesDistintasAlineamiento=which(QueryAlineado!=GenomaAlineado)	## Almacena las posiciones del alineamiento en las que hay diferencias.

	PrimeraPosiciónGenoma=((VentanasPositivas[1]-1)*wsize+SumaResto)

	PosicionesDistintas=PosicionesDistintasAlineamiento+PrimeraPosiciónGenoma-1	## Pequeño ajuste necesario por contarse las mutaciones a partir del 1 y no del 0.

	Mutaciones=QueryAlineado[PosicionesDistintasAlineamiento]	## Almacena las mutaciones.

	Genotipado=data.frame(site=PosicionesDistintas,alt=Mutaciones)	## Se crea un data frame que almacena las posiciones y la mutación en dicha posición.

	CladePos=integer(0)	## Búsqueda en los datos de clados proporcionados. Se almacenan las filas cuyas mutaciones y posiciones coinciden con las del alineamiento.

	Found=integer(0) 							## Se define un vector de iguales caracterí­sticas en el que se almacenarán las posiciones de las mutaciones genotipadas.

	for (n in 1:nrow(Genotipado))		## Bucle donde se comparan las filas del genotipado con las de la base de datos.

		{

			donde = which(CladeData[,2]==Genotipado[n,1] & CladeData[,3]==Genotipado[n,2])	## Almacena las filas que coincidan.

  			if (length(donde>0)) { CladePos = c(CladePos,donde); Found = c(Found,n) }	## Si se encuentran coindicendias, se almacena en CladePos el número de fila de la base de datos
																## donde se ha encontrado. También se almacena un vector con las posiciones en la data frame
																## generada por la función.
		}


	if (length(Found)>0)	## Si hay coincidencias entre el genotipado y la base de datos, se retiran de la data frame de genotipado.

		{

			Genotipado=Genotipado[-Found,]

		}

	if(nrow(Genotipado)>0)	## Si tras retirar las coindicencias aún quedan elementos en la data frame de genotipado, se les asignan clados desconocidos.

		{

			Desc=data.frame(clade=rep("?",nrow(Genotipado)))

			Genotipado=cbind(Desc,Genotipado)

			message("Se han encontrado mutaciones no incluidas en la base de datos introducida. Se recomienda comprobar si está actualizada.
					")

		}

		FoundClades=CladeData[CladePos,]

	

		if(!FullClades)	## Si FullClades equivale a FALSE, el único clado que se mostrará será el más reciente.

			{

				Orden=sort(unique(FoundClades[,"clade"]))

				FoundClades=FoundClades[which(FoundClades[,"clade"]==Orden[length(Orden)]),]

			}

		Resultado=rbind(FoundClades,Genotipado)	## Almacena como resultado la información correspondiente a las mutaciones encontradas.

		if(SiteQuery)					## Si el argumento SiteQuery tiene el valor de TRUE, se actualizan las posiciones de las mutaciones para reflejar la posición en la
									## secuencia genotipada.
			{

			Resultado[,"site"]=(Resultado[,"site"]-PrimeraPosiciónGenoma)+1

			}


		message("Resultados del genotipado:
				")					## Muestra los resultados del genotipado.

		return(Resultado)
	
}