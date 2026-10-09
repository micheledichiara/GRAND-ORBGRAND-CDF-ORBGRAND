// Simulazione di trasmissione su canale BSC, con aggiunta di blocco codificatore e decodificatore GRAND con decisione hard

#define _CRT_SECURE_NO_WARNINGS //dice all'IDE di non segnalare come problema funzioni tipo fopen e fscanf

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define N 128 //lunghezza della sequenza ricevuta e dei pattern di errore
#define NH 64 //numero di righe della matrice parity-check => lunghezza della sindrome

//RAPPRESENTAZIONE SPARSA DELLA MATRICE H
//per ogni riga memorizza gli indici delle 8 colonne contenenti un 1
int H_pos[NH][8];

//CONTATORE DEL NUMERO DI QUERY EFFETTUATE DA GRAND
//ogni query corrisponde a un pattern candidato effettivamente sottoposto al controllo di sindrome
unsigned long long numero_query = 0;

//NUMERO MASSIMO DI QUERY CONSENTITE A GRAND
//se il decoder raggiunge questo numero di query senza aver trovato una codeword valida, la ricerca viene interrotta e il blocco viene classificato come abbandono
#define MAX_QUERY 2500000ULL

//flag che indica se GRAND ha esaurito il numero massimo di query
int abbandono = 0;

//FLAG PER ABILITARE O DISABILITARE LE STAMPE DETTAGLIATE
//1 = stampa completa utilizzata per il debug e la validazione
//0 = output ridotto per le simulazioni su molti blocchi
#define STAMPE_DETTAGLIATE 0

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

//NUOVA VERSIONE SPARSA DEL CALCOLO DELLA SINDROME
//sfrutta il fatto che ogni riga della matrice LDPC H contiene esattamente 8 elementi pari a 1
void calcola_sindrome_sparsa(int Z[], int SZ_sparsa[NH]) {

	int r;
	int j;

	for (r = 0; r < NH; r++) {

		SZ_sparsa[r] = 0;

		//invece di scorrere tutte le N=128 colonne considero soltanto le 8 posizioni in cui la riga di H contiene un 1,
		for (j = 0; j < 8; j++) {

			SZ_sparsa[r] =
				(SZ_sparsa[r] + Z[H_pos[r][j]]) % 2;
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
void prova_pattern(int Y[], int H[NH][N], int Ecan[], int Z[], int SZ[NH], int Xcan[], int Efin[], int* trovato) {
	int i;

	//costruzione della candidata Z = Y ^ Ecan, per togliere al ricevuto il rumore ipotizzato, tramite XOR
	for (i = 0;i < N;i++) {
		Z[i] = Y[i] ^ Ecan[i];
	}

	//CONTROLLO DEL NUMERO MASSIMO DI QUERY

	//prima di effettuare una nuova query verifico se il numero massimo
	//di tentativi consentiti è già stato raggiunto
	if (numero_query >= MAX_QUERY) {

		abbandono = 1; //il decoder deve interrompere la ricerca
		return;
	}

	//ogni volta che una candidata Z viene realmente sottoposta al controllo di sindrome, GRAND ha effettuato una nuova query
	numero_query++;
	
	//calcolo della sindrome della candidata
	//una volta costruita Z, va testata, tramite sindrome
	calcola_sindrome_sparsa(Z, SZ);

	//stampa del tentativo corrente
	if (STAMPE_DETTAGLIATE) {
		stampa_vettore(Ecan, N, "Ecan = ");
		stampa_vettore(Z, N, "Z   = ");
		stampa_vettore(SZ, NH, "SZ   = ");
	}
	//se la sindrome è nulla, la candidata appartiene al codice
	if (sindrome_nulla(SZ) == 1) {
		*trovato = 1; //aggiorno il flag di successo del decoder
		copia_vettore(Xcan, Z, N); //salvo la candidata come codeword finale stimata
		copia_vettore(Efin, Ecan, N); //salvo il pattern di errore che ha portato al successo
	}
}

//funzione più critica del programma, che serve a generare tutti i pattern a seconda del peso
//permette di svincolarsi dal caso rigido di scrivere artigianalmente i pattern per ogni peso
//genera automaticamente i pattern in funzione di N
//uno_da_insierire = variabile che dice quanti 1 mancano da collocare nel pattern; posizione_iniziale = variabile che dice da quale posizione continuare a mettere i prossimi 1, per evitare duplicati 
void genera_pattern(int Y[], int H[NH][N], int Ecan[], int Z[], int SZ[NH], int Xcan[], int Efin[], int* trovato, int uno_da_inserire, int posizione_iniziale) {
	int i; //indice del ciclo di scelta della posizione del prossimo 1

	//la ricorsione deve arrestarsi sia quando viene trovata una codeword valida sia quando viene raggiunto il numero massimo di query
	if (*trovato == 1 || abbandono == 1) {
		return;
	}

	// se il pattern è completo e può essere testato
	if (uno_da_inserire == 0) {
		prova_pattern(Y, H, Ecan, Z, SZ, Xcan, Efin, trovato); //=> il pattern può essere testato
		return; //la chiamata ricorsiva termina
	}

	//CASO RICORSIVO: scelgo dove mettere il prossimo 1
	for (i = posizione_iniziale;i <= N - uno_da_inserire; i++) { //scorre tutte le possibili posizioni in cui mettere il prossimo 1
		//parte da posizione_iniziale perché le precedenti posizioni sono già state decise a livelli ricorsivi superiori
		//arriva a N-uno_da_inserire perché va lasciato abbastanza spazio per gli altri 1 che mancano da collocare, p.e se N=5 e mancano due 1, non si può partire dall'ultima posizione
		Ecan[i] = 1;// mette un uno nell'i-esima posizione

		//si richiama la funzione perché, collocato già un 1, ne manca uno in meno nel processo
		genera_pattern(Y, H, Ecan, Z, SZ, Xcan, Efin, trovato, uno_da_inserire - 1, i + 1); //c'è i+1 perché le prossime scelte devono avvenire nelle sole posizioni successive, evitando duplicati
		Ecan[i] = 0; //annulla la scelta appena fatta, perché, dopo aver cercato tutti i pattern, bisogna tronare indietro e provare un'altra posizione

		//la ricerca deve arrestarsi sia in caso di successo sia in caso di abbandono
		if (*trovato == 1 || abbandono == 1) {
			return;
		}
	}
}

//definizione di una funzione per l'eliminazione gaussiana, prende in ingresso una matrice (A), due indici di riga (r1 e r2), senza restituire nulla (void)
//questa funzione serve a scambiare due righe della matrice
void scambia_righe(int A[NH][N], int r1, int r2) {
	int c; //c serve per scorrere le colonne
	int temp; //temp è una variabile temporanea per lo swap, che serve per non predere un valore durante lo scambio

	//nel ciclo for parto dalla colonna 0, fino a N-1, scorrendo tutte le colonne
	for (c = 0;c < N;c++) {
		temp = A[r1][c]; //salvo temporaneamente l'elemento della riga r1 e colonna c
		A[r1][c] = A[r2][c]; //metto dentro la riga r1 l'elemento che stava nella riga r2
		A[r2][c] = temp; //rimetto nella riga r2 il valore vecchio della riga r1, così le due righe si scambiano
	}
}

//riduzione generale di una matrice binaria tramite eliminazione gaussiana su GF(2)
//la funzione modifica direttamente L, memorizza in pivot_cols le colonne in cui trova i pivot
//e restituisce il numero totale di pivot trovati, cioè il rango della matrice
//È praticamente lo stesso principio già usato nella vecchia eliminazione gaussiana, solo che ora non è più vincolato alla particolare struttura del (5,2).

int riduci_gf2(int L[NH][N], int pivot_cols[NH]) {

	int r_pivot = 0;      //riga nella quale deve essere posizionato il prossimo pivot
	int c_pivot;          //colonna che si sta analizzando per cercare il pivot
	int r;                //indice di riga
	int c;                //indice di colonna
	int r_trovata;        //riga nella quale viene trovato un 1 utilizzabile come pivot

	//scorro tutte le colonne della matrice da sinistra verso destra
	//in forma_sistematica() c'era for (c_pivot = K; c_pivot < N && r_pivot < NH; c_pivot++) perché si cercavano i pivot solo nella parte destra della matrice.
	//qui invece si ha for (c_pivot = 0; c_pivot < N && r_pivot < NH; c_pivot++) quindi inizio dalla prima colonna e analizzo tutte le 128 colonne
	//insomma, invece di cercare l'identità, si cercano le colonne linearmente indipendenti
	for (c_pivot = 0; c_pivot < N && r_pivot < NH; c_pivot++) {

		r_trovata = -1; //all'inizio di ogni nuova colonna
		// -1 è un valore di riferimento che funge da flag per indicare che nello specifico momento non è stata trovata nessuna riga valida

		//cerco un 1 nella colonna corrente, partendo dalla riga pivot corrente
		for (r = r_pivot; r < NH; r++) {

			if (L[r][c_pivot] == 1) { //appena trovata una riga valida, si aggiorna la label r_trovata, se alla fine è ancora -1, significa che quella colonna non può ospitare un nuovo pivot. => non si restituisce errore come nella vecchia funzione
				r_trovata = r;
				break;
			}
		}

		//se nella parte restante della colonna non è presente alcun 1,
		//questa colonna non può essere utilizzata come colonna pivot
		if (r_trovata == -1) {
			continue; //se la colonna non aggiunge un nuovo pivot => si passa alla successiva
		}

		//se il pivot è stato trovato in una riga diversa da quella desiderata,
		//porto la riga trovata nella posizione r_pivot
		if (r_trovata != r_pivot) {
			scambia_righe(L, r_pivot, r_trovata); //richiamo della funzione scambia_righe esattamente come in forma_sistematica()
		}

		//annullo tutti gli altri 1 presenti nella colonna pivot, sia sopra sia sotto il pivot stesso

		for (r = 0; r < NH; r++) { //scorre tutte le righe e non soltanto quelle sotto r_pivot

			if (r != r_pivot && L[r][c_pivot] == 1) {

				for (c = 0; c < N; c++) {
					L[r][c] = (L[r][c] + L[r_pivot][c]) % 2; // Per annullare un altro 1, 1+1=0 (mod-2) => si ottiene la forma desiderata
				}
			}
		}

		//memorizzo quale colonna è stata utilizzata come pivot
		//analogo a una permutazione, => per esempio, se i primi pivot si trovano nelle colonne 1, 4, 5 => pivot_cols[0]=1, pivot_cols[1]=4, pivot_col[3]=5
		//Tutte le colonne non presenti in pivot_cols saranno invece le colonne libere, da cui partità la costruzione di G
		pivot_cols[r_pivot] = c_pivot;

		//il prossimo pivot dovrà essere collocato nella riga successiva
		r_pivot++;
	}

	//il numero di righe pivot trovate coincide con il rango della matrice
	return r_pivot;
}


int main() {

		int x; //singolo bit trasmesso, preso di volta in volta da X[i]
		int n; //rumore
		int y; //singolo bit ricevuto
		int errori = 0;
		double p_flip = 0.004; //probabilità di crossover = 0.1
		double u; //variabile casuale reale nell'intervallo [0,1]
		int a; //serve in seguito per la generazione di rumore
		int X[N]; //parola di codice trasmessa sul canale
		int R[N]; //sequenza di rumore
		int Y[N]; //sequenza ricevuta
		int S[NH]; //sindrome

		int H[NH][N]; //crea la matrice H (64 x 128)
		int L[NH][N]; //copia di lavoro di H, da modificare durante l'eliminazione gaussiana

		int pivot_cols[NH]; //memorizza le colonne in cui vengono individuati i pivot
		int rango; //numero di pivot trovati = rango della matrice H

		int free_cols[N]; //conterrà gli indici delle colonne non pivot
		int num_free; //numero effettivo di colonne libere
		int K; //dimensione del codice

		int r;
		int c;

		int peso; //variabile che rappresenta il peso di Hamming del pattern in corso
		int trovato = 0; //flag che indica sei il decoder ha già trovato una soluzione, all'inzio non ne è stata trovata una

		//definizone dei vettori di lavoro, già visti nella definizione di funzioni, richiamate per il lavoro da compiere
		int Ecan[N];
		int Z[N];
		int SZ[NH];
		int Xcan[N];
		int Efin[N];

	FILE* fp = fopen("H_LDPC_128_wc4_wr8.txt", "r"); //fp è una variabile che rappresenta il file aperto
	//è il collegamento che permette di leggere il file
	//fopen == "open file", è presente "r" perché va aperto in modalità r, ovvero "read" => da non confondere con l'indice del rumore

	//controllo della bontà dell'apertura del file
	//se fp==NULL => il programma non è riuscito ad aprire il file
	if (fp == NULL) {
		printf("Errore nell'apertura del file.\n");
		return 1; // non return 0 al fine di segnalazione di errore
	}

	for (r = 0; r < NH; r++) {
		for (c = 0; c < N; c++) {
			fscanf(fp, "%d", &H[r][c]); //fscanf legge un numero e lo mette nella cella relativa =>H[r][c]
		}
	}

	//dopo la lettura di H, aggiunta di un blocco che copia H in L, sempre tramite cicli annidati
	for (r = 0;r < NH;r++) {
		for (c = 0;c < N;c++) {
			L[r][c] = H[r][c];
		}
	}

	fclose(fp); //chiusura del file

	printf("Matrice H letta da file: \n"); //stampa del titolo
	for (r = 0;r < NH;r++) {
		for (c = 0;c < N;c++) {
			printf("%d ", H[r][c]); //stampa della matrice
		}
		printf("\n");
	}

	//dopo la stampa di H, aggiunta della stampa di L
	printf("Matrice di lavoro L: \n");
	for (r = 0;r < NH; r++) {
		for (c = 0;c < N;c++) {
			printf("%d ", L[r][c]);
		}
		printf("\n");
	}

	//COSTRUZIONE DELLA RAPPRESENTAZIONE SPARSA DI H, serve principalmente per ottimizzare la resa della compilazione dal punto di vista del tempo di elaborazione

	for (r = 0; r < NH; r++) {

		int indice_uno = 0;

		for (c = 0; c < N; c++) {

			if (H[r][c] == 1) {

				H_pos[r][indice_uno] = c;
				indice_uno++;
			}
		}
	}

	printf("Rappresentazione sparsa di H costruita correttamente.\n");
	
	//riduzione della matrice di lavoro L su GF(2)
	//la funzione restituisce il numero di pivot trovati, cioè il rango di H
	rango = riduci_gf2(L, pivot_cols);

	printf("Rango di H calcolato in C: %d\n", rango);

	K = N - rango; //istruzione immediata, implementa il calcolo della lunghezza del codice

	printf("Dimensione del codice K calcolata in C: %d\n", K);

	int is_pivot[N] = { 0 }; //questa inizializzazione porta tutti i 128 elementi a zero.
	//=> all'inizio considero tutte le colonne non pivot
	//Poi marco quelle che riduci_gf2() ha effettivamente trovato


	//Marco con 1 tutte le colonne che riduci_gf2() ha individuato come pivot
	//quando r = 0, prende pivot-cols[0] che contiene l'indice della prima colonna pivot, e si fa is_pivot[pivot_cols[0]]=1
	//Poi per r=1 lo stesso, is_pivot[pivot_cols[1]]=1 eccetera. Il ciclo arriva fino a r<rango => 61 valori
	for (r = 0; r < rango; r++) {
		is_pivot[pivot_cols[r]] = 1;
	}

	//Cerco le colonne che non sono pivot e ne salvo gli indici in free_cols
	//=> operazione opposta
	//Contatore di colonne libere trovate
	num_free = 0;

	for (c = 0; c < N; c++) {

		if (is_pivot[c] == 0) {

			free_cols[num_free] = c;
			num_free++;
		}
	}

	printf("Numero di colonne pivot: %d\n", rango);
	printf("Numero di colonne libere: %d\n", num_free);

	//il numero di colonne libere deve coincidere con la dimensione K del codice
	if (num_free != K) {
		printf("Errore: il numero di colonne libere (%d) non coincide con K (%d).\n", num_free, K);
		return 1;
	}

	printf("Controllo superato: il numero di colonne libere coincide con K.\n");

	//allocazione dinamica della matrice generatrice G
	//G avrà K righe e N colonne, quindi nel caso corrente dimensione 67 x 128
	//G è un puantatore a righe costituite ciascuna da N interi, poiché N=128, ogni riga puntata da G contiene 128 interi 
	int (*G)[N] = calloc(K, sizeof * G); //calloc riserva memoria sufficiente per K righe della stessa dimensionedi una riga di G => verrà allocata una matrice 67x128

	if (G == NULL) {
		printf("Errore nell'allocazione della matrice generatrice G.\n");
		return 1;
	}

	int q;//indice della riga di G, associata a una variabile libera
	int p;//indice utilizzato per scorrere le colonne pivot
	int free_col;//indice fisico della colonna libera corrente

	//costruzione della matrice generatrice G
	//ogni riga di G viene associata ad una delle K colonne libere
	for (q = 0; q < K; q++) {

		//indice fisico della q-esima colonna libera
		free_col = free_cols[q];

		//La variabile libera associata alla riga corrente viene posta uguale a 1
		//le altre colonne libere restano a 0 grazie all'inizializzazione effettuata da calloc
		G[q][free_col] = 1;

		//Determino il valore delle variabili pivot utilizzando le equazioni contenute nella matrice ridotta L

		for (p = 0; p < rango; p++) {

			//Punto fondamentale, una volta scelta la colonna libera associata alla riga corrente di G, leggo nella matrice ridotta L quali valori devono assumere le colonne pivot e li copio nelle corrispondenti posizioni della riga di G.
			//Per costrire ogni riga della matrice generatrice considero una variabile libera alla volta: la pongo a 1, mentre tutte le altre variabili libere rimangono a 0. La matrice ridotta L consente poi di ricavare i valori delle variabili pivot necessari a soddisfare le equazioni di parità. 
			//Per questo, per ogni pivot, il valore da inserire nella riga di G coincide con l'elemento di L posto sulla stessa riga pivot e sulla colonna libera considerata.
			G[q][pivot_cols[p]] = L[p][free_col];
		}
	}

	//CONTROLLO DELLA PARTE DI G CORRISPONDENTE ALLE COLONNE LIBERE

	//durante la costruzione di G, per ogni riga q, è stata posta a 1 solamente la q-esima variabile libera
	//tutte le altre variabili libere sono invece rimaste a 0 grazie all'inizializzazione effettuata con calloc
	//=> se considero solamente le colonne contenute in free_cols, devo quindi ottenere una matrice identità KxK

	int j_controllo; //indice che serve per scorrere le K colonne libere
	int controllo_identita = 1; //Flag di controllo; parto assumendo che il controllo sia corretto
	//se durante la verifica trovo anche un solo elemento errato, questo flag viene posto a 0

	for (q = 0; q < K; q++) {

		for (j_controllo = 0; j_controllo < K; j_controllo++) {

			if (G[q][free_cols[j_controllo]] != (q == j_controllo)) {

				controllo_identita = 0;

				break;
			}
		}

		if (controllo_identita == 0) {
			break;
		}
	}

	if (controllo_identita == 1) {

		printf("Controllo G superato: sulle colonne libere e' presente I_K.\n");
	}
	else {

		printf("Errore nella costruzione di G: il blocco sulle colonne libere non e' I_K.\n");

		free(G);

		return 1;
	}


	//CONTROLLO DI ORTOGONALITA' TRA H E G

	int controllo_ortogonalita = 1;
	int somma;

	for (r = 0; r < NH; r++) {

		for (q = 0; q < K; q++) {

			somma = 0;

			for (c = 0; c < N; c++) {

				somma = (somma + H[r][c] * G[q][c]) % 2;
			}

			if (somma != 0) {

				controllo_ortogonalita = 0;

				printf("Errore di ortogonalita' tra la riga %d di H e la riga %d di G.\n",
					r, q);

				break;
			}
		}

		if (controllo_ortogonalita == 0) {
			break;
		}
	}

	if (controllo_ortogonalita == 1) {

		printf("Controllo di ortogonalita' superato: H * G^T = 0.\n");
	}
	else {

		printf("Errore: la matrice G costruita non soddisfa H * G^T = 0.\n");

		free(G);
		return 1;
	}

	//allocazione dinamica del messaggio informativo M, formato da K bit, => si chiede memoria per K=67 interi
	int* M = malloc(K * sizeof * M);

	if (M == NULL) {
		printf("Errore nell'allocazione del messaggio informativo M.\n");
		free(G);
		return 1;
	}

	//si inizializza il generatore pseudo-casuale
	srand((unsigned int)time(NULL));

	//CONTATORI PER IL TEST SU PIU' BLOCCHI

	int numero_blocchi = 0;       //numero totale di blocchi trasmessi
	int numero_corretti = 0;      //blocchi decodificati correttamente
	int numero_undetection = 0;   //codeword valida trovata, ma diversa da X
	int numero_abbandoni = 0;     //decoder fermato per MAX_QUERY

	//numero di errori totali che si vogliono raccogliere per ogni valore di p
	//la simulazione continua finché la somma tra undetection e abbandoni raggiunge tale valore
	//=> TARGET_ERRORI rappresenta il criterio di arresto della simulazione
	const int TARGET_ERRORI = 100;

	unsigned long long query_totali = 0; //somma delle query effettuate su tutti i blocchi
	unsigned long long query_massime_successo = 0; //massimo numero di query osservato nei blocchi in cui GRAND trova una codeword


	//SIMULAZIONE SU UN NUMERO VARIABILE DI BLOCCHI
	//la simulazione continua finché vengono osservati 100 errori totali, dove per errore totale si intende undetection + abbandono
	//=> il numero di blocchi trasmessi dipende dalle prestazioni di GRAND al valore di p considerato

	while (numero_undetection + numero_abbandoni < TARGET_ERRORI) {

		//ogni iterazione del while corrisponde alla trasmissione
		//e successiva decodifica di un nuovo blocco
		numero_blocchi++;

		printf("\n\n=========================================\n");
		printf("BLOCCO %d\n", numero_blocchi);
		printf("=========================================\n");

		errori = 0;


		//GENERAZIONE CASUALE DI M

		//generazione casuale dei K bit del messaggio informativo
		//ogni posizione può assumere solamente 0 oppure 1
		//=> rand()%2 restituisce direttamente un bit binario
		for (c = 0; c < K; c++) {
			M[c] = rand() % 2;
		}

		//CODIFICA TRAMITE G
		//Per costruire ogni bit della codeword scorro tutti i bit del messaggio informativo e considero il contributo della corrispondente riga di \(G\); la somma viene effettuata modulo 2.
		int i; //indice utilizzato per scorrere i K bit del messaggio informativo

		//CODIFICA DEL MESSAGGIO INFORMATIVO
		//calcolo X = M * G in GF(2)
		for (c = 0; c < N; c++) {

			X[c] = 0; //azzero il bit della codeword prima di calcolarlo

			for (i = 0; i < K; i++) {

				X[c] = (X[c] + M[i] * G[i][c]) % 2;
			}
		}

		//CONTROLLO DELLA CODEWORD X TRAMITE SINDROME
		//Dopo aver codificato M tramite G, X deve essere una parola valida del codice
		//=> Deve soddisfare tutti i vincoli di parità imposti da H
		//=> calcolo la sindrome S_X = X * H^T su GF(2), che deve risultare completamente nulla

		int SX[NH]; //sindrome della codeword X, composta da NH=64 elementi
		int controllo_sindrome_X = 1; //flag; inizialmente assumo che la sindrome sia nulla

		for (r = 0; r < NH; r++) { //ogni riga di H corrisponde ad un vincolo di parità

			SX[r] = 0; //azzero l'elemento corrente della sindrome

			for (c = 0; c < N; c++) { //scorro tutti i 128 bit della codeword

				//calcolo il prodotto scalare mod-2 tra X e la riga r di H
				//=> se il risultato è 0, X soddisfa quello specifico vincolo di parità
				SX[r] = (SX[r] + X[c] * H[r][c]) % 2;
			}

			//se anche un solo elemento della sindrome vale 1,
			//la sequenza X ottenuta dalla codifica non è una parola valida del codice
			if (SX[r] != 0) {

				controllo_sindrome_X = 0;

				printf("Errore: sindrome non nulla alla posizione %d.\n", r);

				break;
			}
		}

		//verifica finale del controllo
		if (controllo_sindrome_X == 1) {

			if (STAMPE_DETTAGLIATE) {
				printf("Controllo codifica superato: la codeword X ha sindrome nulla.\n");
			}
		}
		else {

			printf("Errore nella codifica: X non soddisfa tutti i vincoli di parita' di H.\n");

			free(M);
			free(G);

			return 1;
		}

		//CONTROLLO DELLA SISTEMATICITA' DELLA CODIFICA
		//la matrice G costruita contiene una matrice identità nelle colonne libere
		//=> ogni bit informativo M[q] deve essere riportato direttamente
		//nella posizione della codeword indicata da free_cols[q]

		int controllo_sistematicita = 1; //flag; inizialmente assumo che il controllo sia corretto

		for (q = 0; q < K; q++) {

			//confronto il q-esimo bit informativo con la posizione della codeword
			//associata alla q-esima colonna libera
			if (X[free_cols[q]] != M[q]) {

				controllo_sistematicita = 0;

				printf("Errore di sistematicita': M[%d] diverso da X[%d].\n",
					q, free_cols[q]);

				break;
			}
		}

		if (controllo_sistematicita == 1) {

			if (STAMPE_DETTAGLIATE) {
				printf("Controllo sistematicita' superato: M coincide con X nelle colonne libere.\n");
			}
		}

		else {

			free(M);
			free(G);

			return 1;
		}

		//inserimento nel canale rumoroso con eventuale crossover

		// simulazione bit per bit della trasmissione

		for (i = 0; i < N;i++) {
			x = X[i]; //estrazione del bit trasmesso in posizione i

			//generazione di un valore randomico tra 0 e 1
			// ciò è possibile solo dividendo a per RAND_MAX, cioè estremo dell'intervallo di valori di scelta pseudo-casuale
			a = rand();
			u = (double)a / RAND_MAX;

			//generazione del rumore binario, confrontando r con p
			if (u < p_flip) {
				n = 1;
			}
			else {
				n = 0;
			}

			//calcolo del bit ricevuto tramite XOR tra bit trasmesso e rumore
			y = x ^ n;

			R[i] = n;
			Y[i] = y;


			//aggiornamento del contatore di errori di trasmissione
			if (x != y) {
				errori++;
			}

			//stampa della trasmissione per il bit specifico
			//printf("\ni=%d, x=%d, n=%d, y=%d\n", i, x, n, y);
		}

		//inserimento matrice di controllo di parità e conseguente calcolo della sindrome, per la verifica di presenza di errori nella trasmissione

		for (r = 0;r < NH;r++) {
			S[r] = 0;
			for (c = 0;c < N;c++) {
				S[r] = (S[r] + Y[c] * L[r][c]) % 2;
			}
		}

		//oppure calcola_sindrome(Y, L, S);, si fa anche prima veramente, però lo si lascia per coerenza


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

		printf("\nRumore R: ");
		for (i = 0; i < N; i++) {
			printf(" % d", R[i]);
		}

		printf("\nSequenza ricevuta Y: ");
		for (i = 0; i < N; i++) {
			printf("%d", Y[i]);
		}

		printf("\nErrori totali = %d\n", errori);
		printf("Numero totale di bit trasmessi = %d\n", N);
		printf("Probabilita' di flip p= %f\n", p_flip);

		//inserimento sindrome e verifica del suo valore
		printf("\nSindrome S: ");
		for (r = 0;r < NH;r++) {
			printf("%d", S[r]);
		}


		if (sindrome_nulla(S) == 1) {
			printf("\nLa sequenza ricevuta soddisfa i vincoli di parita' del codice.\n");
		}
		else {
			printf("\nLa sequenza ricevuta non soddisfa i vincoli di parita' del codice.\n");
		}

		trovato = 0;
		numero_query = 0;
		abbandono = 0;

		//ciclo sul peso di Hamming
		//si provano tutti i pattern di errore, partendo dai più probabili, organizzati per peso crescente, fino a N, massimo numero di errori 
		for (peso = 0;peso <= N;peso++) {
			printf("\n=====================\n");
			printf("PATTERN DI PESO %d\n", peso); //stampa del peso corrente di lavorazione
			printf("=====================\n");

			azzera_vettore(Ecan, N); //azzeramento di Ecan, perché bisogna partire da una sequenza pulita
			genera_pattern(Y, H, Ecan, Z, SZ, Xcan, Efin, &trovato, peso, 0); //generazione di tutti i pattern di peso "peso"; c'è 0 = prima posizione da cui è possibile partire c'è peso = quanti 1 vanno collocati, inoltre c'è &trovato per aggiornare la variabile originale, non una copia

			if (trovato == 1 || abbandono == 1) {
				break;
			}
		}

		//BLOCCO FINALE DI STAMPA E VERIFICA DEL DECODER

		printf("\n=====================================\n");
		printf("RISULTATO FINALE DEL DECODER\n");
		printf("=====================================\n");

		if (trovato == 1) {

			if (numero_query > query_massime_successo) {
				query_massime_successo = numero_query;
			}

			//MODIFICA RISPETTO AL VECCHIO PROGETTO:
			//nel progetto toy veniva utilizzata stampa_vettore().
			//Nel nuovo progetto, con N=128, i vettori vengono stampati direttamente
			//disponendo gli elementi su più righe per rendere l'output più leggibile.

			if (STAMPE_DETTAGLIATE) {

				printf("\nCodeword stimata Xcan:\n");
				for (i = 0; i < N; i++) {
					printf("%d ", Xcan[i]);

					if ((i + 1) % 16 == 0) {
						printf("\n");
					}
				}

				printf("\nPattern trovato Ecan finale (Efin):\n");
				for (i = 0; i < N; i++) {
					printf("%d ", Efin[i]);

					if ((i + 1) % 16 == 0) {
						printf("\n");
					}
				}

				printf("\nSindrome finale SZ:\n");
				for (i = 0; i < NH; i++) {
					printf("%d ", SZ[i]);

					if ((i + 1) % 16 == 0) {
						printf("\n");
					}
				}
			}

			printf("\nGRAND ha trovato una parola valida del codice.\n");

			printf("Numero totale di query effettuate: %llu\n", numero_query);

			int uguali = 1;

			for (i = 0; i < N; i++) {

				if (X[i] != Xcan[i]) {
					uguali = 0;
					break;
				}
			}

			if (uguali == 1) {

				numero_corretti++;

				printf("Decodifica corretta: Xcan coincide con la sequenza trasmessa X.\n");
			}
			else {

				numero_undetection++;

				printf("GRAND ha trovato una parola valida del codice, ma diversa dalla sequenza trasmessa X.\n");
			}
		}
		else {

			if (abbandono == 1) {

				numero_abbandoni++;

				printf("\nGRAND ha raggiunto il numero massimo di query.\n");
				printf("Numero totale di query effettuate: %llu\n", numero_query);
				printf("Limite massimo di query: %llu\n", (unsigned long long)MAX_QUERY);
				printf("Decodifica interrotta: il blocco viene classificato come ABBANDONO.\n");
			}

			else {

				printf("GRAND non ha trovato una parola valida tra i pattern testati.\n");
			}
		}

		query_totali += numero_query;
	} // fine ciclo sui blocchi

	
	  //RIEPILOGO DEL TEST SU PIU' BLOCCHI
	
	int numero_errori = numero_undetection + numero_abbandoni;

	//CALCOLO DELLA BLOCK ERROR RATE
	//numero_errori rappresenta il numero complessivo di blocchi errati osservati, ottenuto come somma tra undetection e abbandoni
	//al termine del while, in condizioni normali, numero_errori = TARGET_ERRORI = 100
	//numero_blocchi rappresenta il numero totale di blocchi necessario per osservare tali 100 errori
	//=> BLER = numero_errori / numero_blocchi = 100 / numero_blocchi

	double BLER = (double)numero_errori / numero_blocchi;

	double query_medie =
		(double)query_totali / numero_blocchi;

	printf("\n\n=========================================\n");
	printf("RIEPILOGO SIMULAZIONE\n");
	printf("=========================================\n");

	printf("Probabilita' di crossover p: %.6f\n", p_flip);
	printf("Blocchi totali: %d\n", numero_blocchi);
	printf("Decodifiche corrette: %d\n", numero_corretti);
	printf("Undetection: %d\n", numero_undetection);
	printf("Abbandoni: %d\n", numero_abbandoni);
	printf("Errori totali: %d\n", numero_errori);

	printf("BLER: %.8f\n", BLER);

	printf("Query totali: %llu\n", query_totali);
	printf("Query medie per blocco: %.2f\n", query_medie);
	printf("Massimo numero di query tra i successi: %llu\n", query_massime_successo);

	free(M);
	free(G);
	
	return 0;

}