// Simulazione di trasmissione BPSK su canale AWGN, con blocco codificatore e successiva decodifica ORBGRAND

#define _CRT_SECURE_NO_WARNINGS //dice all'IDE di non segnalare come problema funzioni tipo fopen e fscanf

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define N 5 //lunghezza della sequenza ricevuta e dei pattern di errore
#define K 2 //lunghezza del messaggio informativo
#define NH 3 //numero di righe della matrice parity-check => lunghezza della sindrome

//definizione di una funzione che serve a stampare un qualsiasi vettore binario
//invece di riscrivere ogni volta lo stesso ciclo for per ogni vettore, creo una funzione generale
//V[] = vettore da stampare; dim = la sua dimensione, titolo[] = testo da stampare prima del vettore
void stampa_vettore(int V[], int dim, const char titolo[]) { //la funzione non restituisce un valore
	int i; //indice del ciclo, per scorrere le posizioni del vettore da 0 a dim-1

	printf("%s", titolo); //stampa del testo introduttivo
	for (i = 0;i < dim;i++) { //ciclo da 0 a dim-1, per stampare tutti gli elementi del vettore
		printf("%d", V[i]); //stampa dell'elemento del vettore in corso
	}
	printf("\n");
}

//definizione di una funzione che serve a portare una sequenza nella forma 000...0
//utile in corso d'opera per la sequenza Ecan, perchéogni ricerca deve partire da un pattern pulito
//V[] = vettore da azzerare; dim = la sua dimensione
void azzera_vettore(int V[], int dim) {
	int i;

	for (i = 0;i < dim;i++) { //il ciclo scorre tutte le posizioni del vettore e mette ogni elemento a zero, altrimenti Ecan potrebbe contenere valori dei tentativi precedenti
		V[i] = 0;
	}
}

//definizione di una funzione che non restituisce nulla, che serve a copiare il contenuto di un vettore in un altro
//utile se si trova una soluzione valida e si vuole salvare la code word finale in Xcan e il pattern finale in Efin
//destinazione[] = vettore dove salvare i valori; dorgente[] = vettoreda cui si prendono, dim  = numero di elementi da copiare
void copia_vettore(int destinazione[], int sorgente[], int dim) { //per ogni posizione si prende l'elemento sorgente e si salva in destinazione
	int i;

	for (i = 0;i < dim;i++) {
		destinazione[i] = sorgente[i];
	}
}

//definizione di una funzione che realizza il calcolo della sindrome SZ = ZH^T
//Z[] = parola candidata; H[NH][N] = matrice di controllo di parità(nel codice è L, nell'effettivo la sua copia di lavoro, nel codice più importante)
//SZ[NH] = vettore sindrome da riempire
void calcola_sindrome(int Z[], int H[NH][N], int SZ[NH]) {
	int r; //indice di riga
	int c; //indice di colonna

	//ogni riga di H rappresenta un vincolo di parità, quindi ogni riga produce una componente della sindrome
	for (r = 0;r < NH;r++) {
		SZ[r] = 0; //prima di accumulare i contributi della riga r-esima, si azzera a componente della sindrome, perché altrimenti dentro SZ[r] ci possono essere residui di calcolo precedente
		for (c = 0;c < N;c++) {
			SZ[r] = (SZ[r] + Z[c] * H[r][c]) % 2; //prodotto tra il simbolo di Z in posizione c per il corrispondente simbolo della riga r della matrice
		}
	}
}

//definizione di una funzione che serve a controllare se la sindrome è nulla, restituisce 0 o 1, in base al responso (1=sì, 0=no, logica booleana)
int sindrome_nulla(int SZ[NH]) {
	int r; //indice che scorre le componenti del vettore sindrome

	//scorrendo le componenti della sindrome, appena trova anche solo una componenente != 0, so già che la sindrome non è nulla => return 0
	for (r = 0;r < NH;r++) {
		if (SZ[r] != 0) {
			return 0;
		}
	}

	return 1; //se si esce dal ciclo senza aver trovato componenti non nulle, allora tutte sono 0
} //else superfluo, perché basta che il ciclo finisca senza trovare alcun 1 affinché la sindrome sia nulla e non servono condizioni

//definizione  di una funzione che prende un pattern candidato già costruito e lo mette alla prova vedendo se sia quello giusto
//Ecan[] = pattern di errore candidato; Xcan[] = code word finale salvata dal decoder
//Efin[] = pattern finale salvato dal decoder; *trovato = puntatore alla variabile che indica se il decoder ha già trovato una soluzione
void prova_pattern(int Yhard[], int H[NH][N], int Ecan[], int Z[], int SZ[NH], int Xcan[], int Efin[], int* trovato) {
	int i;

	//costruzione della candidata Z = Y ^ Ecan, per togliere al ricevuto il rumore ipotizzato, tramite XOR
	for (i = 0;i < N;i++) {
		Z[i] = Yhard[i] ^ Ecan[i];
	}

	//calcolo della sindrome della candidata
	//una volta costruita Z, va testata, tramite sindrome
	calcola_sindrome(Z, H, SZ);

	//stampa del tentativo corrente
	stampa_vettore(Ecan, N, "Ecan = ");
	stampa_vettore(Z, N, "Z   = ");
	stampa_vettore(SZ, NH, "SZ   = ");

	//se la sindrome è nulla, la candidata appartiene al codice
	if (sindrome_nulla(SZ) == 1) {
		*trovato = 1; //aggiorno il flag di successo del decoder
		copia_vettore(Xcan, Z, N); //salvo la candidata come codeword finale stimata
		copia_vettore(Efin, Ecan, N); //salvo il pattern di errore che ha portato al successo
	}
}

//Definisco una funzione che rstituisce un valore reale. In questa fase voglio ottenere una normale standard, => media nulla e varianza 1.
double genera_gaussiana_standard(void) {

	const double PI = 3.14159265358979323846;


	double U1; //generazione degli ingressi della funzione
	double U2;
	double Z;

	//generazione di due variabili uniformi nell'intervallo (0,1)

	U1 = ((double)rand() + 1.0) / ((double)RAND_MAX + 2.0); //la presenza di +1.0 e +2.0 sono importanti per ottenere calori strettamente compresi tra 0 e 1
	U2 = ((double)rand() + 1.0) / ((double)RAND_MAX + 2.0);

	//in particolare è importante per U1, infatti ln(U1) fa si che U1 debba essere !=0

	//trasformazione di Box-Muller, che trasforma due variabili uniformi in una normale

	Z = sqrt(-2.0 * log(U1)) * cos(2.0 * PI * U2);

	return Z; //per restituire il campione al main
}

//definizione di una funzione per ordinare le affidabilità, universale per poter essere usata anche a progetto completo
void ordina_affidabilita(double Affidabilita[], int pi[]) { //riceve i due vettori Affidabilita[] e pi[], è void dunque modifica i valori di pi[]


	// Affidabilita[] contiene i moduli degli LLR e verrà soltsnto letto, mentre pi[] gli indici delle posizioni fisiche e verrà modificato
	int i; // Posizione di pi che va sistemata
	int j; //serve per cercare, nelle posizioni successive, l'affidabilità più piccola
	int indice_min; //Memorizza la posizione dentro pi nella quale, fino a quel momento, si ha trovato l'affidabilità minima
	int temp; //per scambiare elementi

	for (i = 0;i < N - 1;i++) { //=> selection sort, va decisio quale indice mettere in pi[i], si cerca tra p[i], pi[i+1] ecc..., quello associato all'affidabilità più piccola

		// si parte da i=0, => va trovato il bit meno affidabile di tutti e metterne l'indice in pi[0]
		//poi per i=1 si cercherà il meno affidabile tra quelli restanti e così via
		// //c'è N-1 perché l'ultimo bit rimasto non ha bisogno di cercare posizione

		indice_min = i; //suppongo che il minimo sia quello che si trova nella poszione i di pi. => poi controllo gli altri e vedo se ne trovo uno più piccolo

		//si controllano tutti gli elementi successivi
		//Dato che si ha già assunto quello in poszione i come minimo provvisorio, si parte da J=i+1
		// con i=0 si ha j=1,2,3,4 => va confrontato l'elemento attuale con ciò che segue
		for (j = i + 1;j < N;j++) {
			if (Affidabilita[pi[j]] < Affidabilita[pi[indice_min]]) { // non si confronta l'affidabilità, perché pi potrebbe essere già stato riordinato.

				indice_min = j; //si seguono gli indici contenuti in pi
				//P.e pi = 0 1 2 3 4 ,  per j=1 si ha pi[j]=pi[1]=1; => Affidabilita[pi[j]] diventa Affidabilita[1]
				//Dall'altra parte indice_min=0 => pi[indice_min]=pi[0]=0 e Affidabilita[pi[indice_min]] diventa Affidabilita[0]
				// nell'esempio si aveva Affidabilità[pi[j]] => Affidabilita[1] = 1.90 e Affidabilita[pi[indice_min]] => Affidabilita[0]=4.60
				// => Chiedo 1.90 < 4.60 ? sì, allora eseguo indice_min=j => indice_min =1 => il minimo trovato non è più in poszione 0 di pi, ma 1
				//itero il processo
			}
		}

		//qui chiedo se il valore minimo trovato era già nella posizione corretta o va spostato
		//P.e se si aveva i=0 e poi indice_min=4 => 4 != 0 => va spostato, lo scambio si fa con la variabile temp
		if (indice_min != i) {

			temp = pi[i];
			pi[i] = pi[indice_min];
			pi[indice_min] = temp;
		}
	}
	//qui si usa una semplice selection sort, ma non si confronta Affidabilità[j] direttamente
	//si confronta Affidabilità[pi[j]], perché si riordinano gli indici, non il vettore di affidabilità
	//in altre parole non si vuole trovare il valore dei moduli * fattore di scala, ma degli interi ordinati, che corrispondono alle posizioni fisiche dei bit
	//occorre, d'altra parte, fare attenzione agli indici in C; infatti il C utilizza lo 0, non si avrà p.e 1,2,3,4,5 ; ma 0,1,2,3,4 .
	// => nella stampa finale si aggiungerà pi[i] + 1;
}


//Definizione di una funzione per trasformare il pattern di rango nel pattern fisico
void mappa_pattern(int Erank[], int Ecan[], int pi[]) { //Riceve il pattern dei ranghi e quello fisico da costruire insieme alla permutazione \pi, non restituendo nulla perché serve a modificare Ecan

	int k; //indice tipico del rango, k=0=>rango 1, k=1 => rango 2 [...]

	azzera_vettore(Ecan, N); //prima di costruire un nuovo pattern voglio Ecan=00000, altrimenti si rischia di conservare gli uni appartenenti al pattern precedente

	for (k = 0;k < N;k++) { //ciclo for per scorrere i ranghi

		if (Erank[k] == 1) { //=> P.e se si ha Erank=11000, l'if sarà vero per k=0 e k=1

			Ecan[pi[k]] = 1; //trasformazione da rango a posizione fisica
			//P.e con pi= 3 2 1 4 0 => in realtà (4,3,2,5,1), per k=0 => pi[0]=3 => Ecan[3]=1=> mette 1 in posizione 3 (in C, => posizione 4 matematica)
			//Poi K=1 => pi[1]=2 => Ecan[2]=1 => Ecan = 0 0 1 1 0 
		}
	}
}

void genera_pattern_logistici(int peso_rimanente, int rango_iniziale, int Erank[], int Ecan[], int pi[], int Yhard[],
	int H[NH][N], int Z[], int SZ[NH], int Xcan[], int Efin[], int* trovato) {

	

	//peso_rimanente è una variabile che indica lo stato della generazione, p.e. se all'inizio si ha peso_rimanente=5, scelto il rango 1 si avrà 5-1=4, poi, scelto il rango 4 => peso_rimanente==0, a quel punto vi sarà una combinazione valida
	//rango_iniziale --> serve a impedire di utilizzare due volte lo stesso rango ed evitare duplicati
	int rango;

	if (*trovato == 1) {
		return;
	}
	
	if (peso_rimanente == 0) { //il logistic Weight sarà stato costruito completamente, Erank contiene quindi un pattern valido nel dominio dei ranghi, => lo trasformo nelle corrispondenti posizioni fisiche tramite \pi

		mappa_pattern(Erank, Ecan, pi); //=> attuo la trasformazione tramite \pi e stampo il pattern trovato con la funzione stampa_vettore

		stampa_vettore(Erank, N, "Erank = ");

		prova_pattern(Yhard, H, Ecan, Z, SZ, Xcan, Efin, trovato);

		printf("\n");

		return;
	}

	for (rango = rango_iniziale; rango <= N; rango++) { //Ciclo for con cui attuo uno alla volta tutti i ranghi ancora utilizzabili
		
		//se scelgo un rango, lo inserisco nel pattern, sottraggo il suo valore dal Logistic Weight ancora da costruire, provo ricorsivamente a completare il pattern utilizzando soltanto ranghi successivi e, terminato quel tentativo, annullo la scelta per provarne un'altra
		//qui rango rappresenta il rango di affidabilità che stiamo valutando come possibile elemento del pattern. Non parte sempre da 1, ma da rango_iniziale
		//Un rango può essere selezionato solo se il suo contributo al Logistic Weight non supera il peso ancora da costruire. Se, p.e., peso_rimanente=3, non avrebbe senso scegliere rango=4, perché la somma supererebbe già il Logistic Weight desiderato.

		if (rango <= peso_rimanente) {

			Erank[rango - 1] = 1; //Seleziono provvisoriamente il rango corrente inserendo un 1 nella corrispondente posizione di Erank. Poiché i ranghi matematici vanno da 1 a N mentre gli indici C vanno da 0 a N-1, il rango r viene memorizzato nella posizione r-1.

			//Una volta scelto il rango corrente, richiamo ricorsivamente la funzione per completare il pattern. Il peso ancora da costruire viene ridotto del valore del rango scelto: nuovo_peso = peso_rimanente - rango.
			//La ricerca successiva parte da rango + 1, così lo stesso rango non può essere utilizzato nuovamente e le combinazioni vengono costruite in ordine strettamente crescente, evitando duplicati.
			//Erank mantiene memoria dei ranghi già selezionati durante il percorso ricorsivo; Ecan e pi serviranno per la successiva conversione del
				//pattern di rango nelle effettive posizioni fisiche.

			genera_pattern_logistici(peso_rimanente - rango, rango + 1, Erank, Ecan, pi, Yhard, H, Z, SZ, Xcan, Efin, trovato);

			Erank[rango - 1] = 0; //Terminata l'esplorazione ricorsiva di tutti i pattern che contengono il rango corrente, annullo la scelta riportando a 0 la corrispondente posizione di Erank
			//Questo è il passo di backtracking: consente di ritornare allo stato precedente del pattern e provare, nell'iterazione successiva del for,
			//un rango differente senza conservare gli 1 appartenenti al ramo appena esplorato

			if (*trovato == 1) {
				return;
			} //la stessa funzione in "AWGN" si trasforma da generatore di prova a generatore del decoder vero

		}
	}
}


int main() {
	int i = 0;
	int j = 0;
	int X[N]; //parola di codice trasmessa sul canale

	//definisco tre vettori, X_i è ogni bit codificato
	double Sbpsk[N];
	double RumoreCanale[N];
	double Ysoft[N];

	int M[K]; //definizione messaggio informativo
	int G[K][N]; //definizione matrice generatrice
	int H[NH][N]; //matrice di controllo di parità
	int S[NH]; //sindrome
	int L[NH][N]; //aggiunta della matrice L, copia effettiva della matrice H da file, per la sua manipolazione
	int PT[NH][K];//dichiarazione delle nuove matrici (P TRASPOSTA, P, MATRICE IDENTITA')
	int P[K][NH];
	int rr; //indice di riga
	int cc; //indice di colonna
	int trovato = 0; //flag che indica sei il decoder ha già trovato una soluzione, all'inzio non ne è stata trovata una

	int wL; //Logistic Weight corrente

	int wL_max = N * (N + 1) / 2; //massimo Logistic Weight possibile
	//queste due variabili servono perché il decoder non proverà più solatanto \omega_L=3, ma si generalizza automaticamente \omega_L=0,1,...,\omega_L_max

	int Yhard[N]; //Vettore della decisione hard
	double LLR[N]; //Log-Likelihood Ratio
	double Affidabilita[N]; //Modulo degli LLR
	int pi[N]; //Permutazione degli indici secondo affidabilità

	//definizone dei vettori di lavoro, già visti nella definizione di funzioni, richiamate per il lavoro da compiere
	
	int Erank[N]; //Rappresenta il pattern nel dominio dei ranghi
	
	//Per esempio, Erank=11000 significa "ribalta rango 1 e rango 2"
	
	int Ecan[N];
	int Z[N];
	double Zgauss;     // campione gaussiano standard N(0,1)
	int SZ[NH];
	int Xcan[N];
	int Efin[N];

	//Messaggio informativo
	M[0] = 1;
	M[1] = 1;

	//Parametri del codice
	double Rc = (double)K / N; // il (double) è essenziale perché N e K sono definiti come costanti intere, dunque 2/5 verrebbe 0 se non si entra nel dominio double

	//energia del simbolo BPSK normalizzata
	double Es = 1.0;

	//Rapporto Eb/N0 scelto per la simulazione
	double EbN0_dB = 4.0; //si aggiunge Eb/N0 che conduce a \sigma per analizzare il rumore del canale
	//tale per cui la variabile Z distribuita con PDF normale standard sarà esattamente un elemento con cui si definisce il rumore additivo, con n=\sigma*Z
	//fissato Eb/N0 a 4 dB come valore di prova
	// con Eb/N0=4dB, R=0.4 ed Es=1, => \sigma=0.705, come previsto da calcoli carta e penna
	//con un tale \sigma si ha un rumore abbastanza consistente per essere osservabile

	//conversione da dB a unità assoluta
	double gamma_b = pow(10.0, EbN0_dB / 10.0);

	//varianza teorica del rumore AWGN
	double sigma_quadro = Es / (2.0 * Rc * gamma_b);

	//deviazione standard del rumore
	double sigma = sqrt(sigma_quadro);


	FILE* fp = fopen("lettura_h_nosist.txt", "r"); //fp è una variabile che rappresenta il file aperto
	//è il collegamento che permette di leggere il file
	//fopen == "open file", è presente "r" perché va aperto in modalità r, ovvero "read" => da non confondere con l'indice del rumore

	//controllo della bontà dell'apertura del file
	//se fp==NULL => il programma non è riuscito ad aprire il file
	if (fp == NULL) {
		printf("Errore nell'apertura del file.\n");
		return 1; // non return 0 al fine di segnalazione di errore
	}

	for (rr = 0; rr < NH; rr++) {
		for (cc = 0; cc < N; cc++) {
			fscanf(fp, "%d", &H[rr][cc]); //fscanf legge un numero e lo mette nella cella relativa =>H[rr][cc]
		}
	}

	//dopo la lettura di H, aggiunta di un blocco che copia H in L, sempre tramite cicli annidati
	for (rr = 0;rr < NH;rr++) {
		for (cc = 0;cc < N;cc++) {
			L[rr][cc] = H[rr][cc];
		}
	}

	fclose(fp); //chiusura del file

	printf("Matrice H letta da file: \n"); //stampa del titolo
	for (rr = 0;rr < NH;rr++) {
		for (cc = 0;cc < N;cc++) {
			printf("%d ", H[rr][cc]); //stampa della matrice
		}
		printf("\n");
	}

	//dopo la stampa di H, aggiunta della stampa di L
	printf("Matrice di lavoro L: \n");
	for (rr = 0;rr < NH; rr++) {
		for (cc = 0;cc < N;cc++) {
			printf("%d ", L[rr][cc]);
		}
		printf("\n");
	}

	//ELIMINAZIONE GAUSSIANA
	//per una tale matrice, risulta semplice e manuale la sua implementazione, senza utilizzare specifiche librerie in C
	//è sufficiente applicare alla matrice R1 <- R1+R2 (mod-2)
	//passando da 1 0 1 1 0 a  1 1 1 0 0

	printf("\nEliminazione di Gauss binaria: R1 = R1+ R2 (mod-2)\n");

	//tramite ciclo for, si agisce sulle colonne della matrice, e nell'assegnazione si coinvolge la sola prima riga, => r=0, nell'assegnazione, la cui definizione presuppone
	//il contributo della seconda riga (R1 <- R1+R2), dunque r=1
	//si sommano gli elementi delle due righe e si conserva il risultato nella prima riga => l'implementazione in C presuppone l'assegnazione (=)

	for (cc = 0;cc < N;cc++) {
		L[0][cc] = (L[0][cc] + L[1][cc]) % 2;
	}

	printf("\nMatrice L dopo la riduzione:\n");
	for (rr = 0;rr < NH;rr++) {
		for (cc = 0;cc < N;cc++) {
			printf("%d ", L[rr][cc]);
		}
		printf("\n");
	}

	//estrazione di PT
	//dato che L è già nella forma |P^T|I(n-k)|, PT sta nelle prime 2 colonne, quindi serve un ciclio for per c<K
	for (rr = 0; rr < NH;rr++) {
		for (cc = 0;cc < K;cc++) {
			PT[rr][cc] = L[rr][cc];
		}
	}

	//estrazione di P tramite trasposizione di PT
	//ciò che serve è rendere le righe di PT come colonne di P, quindi nell'associazione (=), l'ordine di c e r cambia membro a membro
	for (rr = 0;rr < K;rr++) {
		for (cc = 0;cc < NH;cc++) {
			P[rr][cc] = PT[cc][rr];
		}
	}

	//volendo ottenere G, di dimensione 2*5, stando alla sua definizione, occorre azzerare la matrice
	//per poi, in base al blocco in cui ci si trova, tradotto c<COLONNE V r<K V c<RIGHE, G è uguale a una matrice identità aggiancata a P
	//lavorando su G, si hanno 2 righe e 5 colonne => K righe e COLONNE colonne
	//nei cicli for r<K sempre, lavorando in G, effettivamente la definizione di G si basa su composizione di I e P in colonne
	// => occorre lavorare per diversi valori dell'indice c

	//azzeramento matrice G, c<COLONNE perché coinvolge tutta la matrice
	for (rr = 0;rr < K;rr++) {
		for (cc = 0;cc < N;cc++) {
			G[rr][cc] = 0;
		}
	}


	//mettere la matrice identità
	//è presente G[r][r], questo perché gli uni si trovano nella diagonale principale, laddove indice di riga e colonna sono uguali
	//=>mi basta lavorare col solo indice riga
	for (rr = 0;rr < K;rr++) {
		G[rr][rr] = 1;
	}

	//copiare P dentro G
	//c<RIGHE, essendo RIGHE=3, c (indice di colonna) deve avere massimo valore in 3, cioè RIGHE
	//aldilà della confusione che può suscitare si mantiene ""RIGHE"" per coerenza col resto del programma
	for (rr = 0;rr < K;rr++) {
		for (cc = 0; cc < NH; cc++) {
			G[rr][K + cc] = P[rr][cc]; //presente G[r][K+c] perché si copia P dopo il blocco identità
		}
	}

	//stampa di PT, P e G

	printf("\nMatrice PT:\n");
	for (rr = 0;rr < NH;rr++) {
		for (cc = 0;cc < K;cc++) {
			printf("%d ", PT[rr][cc]);
		}
		printf("\n");
	}

	printf("\nMatrice P:\n");
	for (rr = 0;rr < K;rr++) {
		for (cc = 0;cc < NH;cc++) {
			printf("%d ", P[rr][cc]);
		}
		printf("\n");
	}

	printf("\nMatrice generatrice G:\n");
	for (rr = 0;rr < K;rr++) {
		for (cc = 0;cc < N;cc++) {
			printf("%d ", G[rr][cc]);
		}
		printf("\n");
	}


	//blocco codificatore

	for (j = 0;j < N;j++) {
		X[j] = 0;
		for (i = 0;i < K;i++) {
			X[j] = (X[j] + M[i] * G[i][j]) % 2;
		}
	}

	//CANALE AWGN 

	srand((unsigned int)time(NULL));

	//dopo la verifica sui 10000 campioni, scelta la convenzione 0 -> +sqrt(Es); 1 -> -sqrt(Es),
	//applico un nuovo ciclo con cui, a partire dai bit, ottengo un set di valori reali (i simboli BPSK) appartenenti a {-1,+1}
	for (i = 0;i < N;i++) {

		//modulazione BPSK
		//si fa tramite if, essendo che a un valore ne faccio associare un altro in dipendenza dal fatto che sia 0 o 1
		if (X[i] == 0) {
			Sbpsk[i] = sqrt(Es); //propriamente la convenzione prevede 0 -> sqrt(Es) , 1 -> -sqrt(Es)
		}
		else {
			Sbpsk[i] = -sqrt(Es);
		}

		//aggiunta del rumore AWGN
		Zgauss = genera_gaussiana_standard();

		//scalatura secondo la \sigma del canale
		RumoreCanale[i] = sigma * Zgauss; //si fa ciò che è stato verificato, ossia Zi con PDF N(0,1) e n_i= \sigma * Z_i, quindi n_i con PDF N(0,\sigma^2)

		//=> per ciascuno dei cinque simboli viene generato un nuovo campione indipendente di rumore 

		//costruisco la Ysoft
		Ysoft[i] = Sbpsk[i] + RumoreCanale[i];

		//dopo aver costruito il singolo valore Ysoft[i], applico la condizione:
		//y_i >= 0 => 0, y_i < 0 => 1
		//=> hard decision => Yhard[i]
		if (Ysoft[i] >= 0) {
			Yhard[i] = 0;
		}
		else {
			Yhard[i] = 1;
		}

		//calcolo degli LLR, viene computato attraverso la relazione:
		// L_i = [ 2 * sqrt(Es) * y_i ] / (\sigma)^2

		LLR[i] = (2.0 * sqrt(Es) * Ysoft[i]) / sigma_quadro;

		//per ricavare l'affidabilità del bit considero il modulo dell'LLR, tramite fabs(), funzione che restituisce il valore assoluto di un double
		Affidabilita[i] = fabs(LLR[i]);

		//nel caso di modulazione standard con Es unitario, la relazione si riduce a LLR = 2 volte Ysoft_i / (\sigma quadro)

		//inizializzazione della permutazione
		pi[i] = i;

	}

	//ordino dal bit meno affidabile al più affidabile tramite la funzione ordina_affidabilita()
	ordina_affidabilita(Affidabilita, pi);




	//stampa delle sequenze M,X,N,Y, del numero totale di errori e probabilità
	//necessario anche qui ciclo for per definire la sequenza posizione per posizione, in maniera corretta
	printf("\nMessaggio informativo M: ");
	for (i = 0;i < K;i++) {
		printf("%d", M[i]);
	}

	printf("\nSequenza codificata X: ");
	for (i = 0;i < N;i++) {
		printf("%d", X[i]);
	}

	printf("\n===============================================\n");
	printf("TRASMISSIONE BPSK SU CANALE AWGN\n");
	printf("===============================================\n");

	printf("X binaria: ");
	for (i = 0; i < N; i++) {
		printf("%d ", X[i]);
	}

	printf("\nS BPSK:    ");
	for (i = 0; i < N; i++) {
		printf("%+.4f ", Sbpsk[i]);
	}

	printf("\nRumore:    ");
	for (i = 0; i < N; i++) {
		printf("%+.4f ", RumoreCanale[i]);
	}

	printf("\nY soft:    ");
	for (i = 0; i < N; i++) {
		printf("%+.4f ", Ysoft[i]);
	}

	printf("\n===============================================\n");
	printf("HARD DECISION E AFFIDABILITA'\n");
	printf("===============================================\n");

	printf("Y hard: ");
	for (i = 0; i < N; i++) {
		printf("%d ", Yhard[i]);
	}

	printf("\nLLR:    ");
	for (i = 0; i < N; i++) {
		printf("%+.4f ", LLR[i]);
	}

	printf("\n|LLR|:  ");
	for (i = 0; i < N; i++) {
		printf("%.4f ", Affidabilita[i]);
	}

	printf("\n\nPermutazione pi (dal meno affidabile al piu' affidabile): ");

	for (i = 0; i < N; i++) {
		printf("%d ", pi[i] + 1);
	}

	printf("\n");



	trovato = 0;


	//ORBGRAND esplora automaticamente i Logistic Weight in ordine crescente
	for (wL = 0; wL <= wL_max; wL++) {

		printf("\n=====================================\n");
		printf("LOGISTIC WEIGHT wL = %d\n", wL);
		printf("=====================================\n");

		//ogni nuovo Logistic Weight deve partire da un pattern di rango nullo
		azzera_vettore(Erank, N);

		//generazione e prova dei pattern aventi Logistic Weight wL
		genera_pattern_logistici(wL,1,Erank,Ecan,pi,Yhard,L,Z,SZ,Xcan,Efin, &trovato);
		//questa chiamata, ad ogni iterazione, presenta: peso_rimanente = wL e rango_iniziale=1
		//=> per wL=0 viene provato il pattern nullo; se fallisce, il for passa a wL=1, poi 2 eccetera

		//appena viene trovata una codeword valida il decoder si arresta
		if (trovato == 1) {
			break;
		}
	}

	//blocco finale di stampa
	printf("\n=====================================\n");
	printf("RISULTATO FINALE DEL DECODER\n");
	printf("=====================================\n");

	if (trovato == 1) {
		stampa_vettore(Xcan, N, "Codeword stimata Xcan: "); //candidata finale
		stampa_vettore(Efin, N, "Pattern trovato Ecan finale (Efin): "); //pattern di errore che ha reso valida la candidata
		stampa_vettore(SZ, NH, "Sindrome finale SZ: "); //relativa sidnrome, nulla
		printf("GRAND ha trovato una parola valida del codice.\n");

		int uguali = 1;
		for (i = 0;i < N;i++) {
			if (X[i] != Xcan[i]) {
				uguali = 0;
				break;
			}
		}

		if (uguali == 1) {
			printf("Decodifica corretta: Xcan coincide con la sequenza trasmessa X.\n");
		}
		else {
			printf("GRAND ha trovato una parola valida del codice, ma diversa dalla sequenza trasmessa X.\n");
		}
	}
	else {
		printf("GRAND non ha trovato una parola valida tra i pattern testati.\n");
	}



	return 0;

}