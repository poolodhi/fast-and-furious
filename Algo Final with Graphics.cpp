// A C Program to demonstrate adjacency list representation of graphs

#include <stdio.h>
#include <stdlib.h>
//#include<iostream.h>
#include<iostream>
#include<limits.h>
#include<math.h>
#include<graphics.h>
//#include<winbgim.h>
//#include<dos.h>
#include<conio.h>
#include<string.h>

#define INF INT_MAX
#define dist_factor 1000
#define Walk_Speed 2
#define Car_Speed 16.2
#define Public_Speed 10.8

#define MAX 1000
#define vertex 112		//plz keep the value one greater than number of nodes eg:-mine were 14 nodes
#define range 10000

#define SCALE 1.5  //Shift from 640*480 to 960*720


//int dist[MAX];//shoretest distace of each vertex in cw path
int V = vertex,k=0,cnode=0,pnode=0;
int visit[vertex],pathnodes[MAX];
int adjnodes[MAX];
int sizeadjnodes;
//int pred[MAX];


using namespace std;



struct Neigh_Node
{
    int dest;
    struct Neigh_Node* next;
    char *Road_Name;
    int distance;       //Road Distance
    int w_time;         //Walking Time
    int c_time;         //Car Traveling Time
    int p_time;         //Public Transport Traveling Time
    int n_xcord;
    int n_ycord;
    int t_density;
};

struct Vertex_List
{
    struct Neigh_Node *head;
    char *Place_Name;
    char Place_Type;
    int v_xcord;
    int v_ycord;
};

struct Graph
{
    int V;
    struct Vertex_List* array;
};
int c=0;

// A utility function to create a new Neighbour Node node
struct Neigh_Node* newNeigh_Node(int value, char *rname, int dist, int wt, int ct, int pt, int t_den, int x, int y)
{
    struct Neigh_Node* node = new Neigh_Node;
    node->dest = value;
    node->next = NULL;
    node->Road_Name = new char[strlen(rname)+1];
    strcpy( node->Road_Name, rname);
    node->distance = dist;
    node->w_time = wt;
    node->c_time = ct;
    node->p_time = pt;
    node->n_xcord = x;
    node->n_ycord = y;
    node->t_density=t_den;
    return node;
}

void addVertex(Vertex_List *a, int index, char const *name, char type, int x, int y)
{
	(a+index)->head = NULL;
    (a+index)->Place_Name = new char[strlen(name)+1];
    strcpy( (a+index)->Place_Name , name );
    (a+index)->Place_Type = type;
    (a+index)->v_xcord = x;
    (a+index)->v_ycord = y;
}

void init_Vertex_List(Vertex_List *a);
void Add_Roads(Graph* graph);


void Draw_Graph(Graph* graph);
void Place_Vertex(int x, int y, char type);
void Draw_Road(int x1, int y1, int x2, int y2, char type,int color);
void Draw_Graph(Graph* graph, int path[], int path_len, int color);
void Display_Name(int x, int y, char* string, int color);

struct Graph* createGraph(int V)
{
    struct Graph* graph = new Graph;
    graph->V = V;


    graph->array = new Vertex_List[V+1];

    init_Vertex_List(graph->array);

    return graph;
}


// Adds an edge to an undirected graph
void addEdge(struct Graph* graph, int src, int dest, char *rname, int dist, char type, int t_den)
{
    // Add an edge from src to dest.  A new node is added to the adjacency
    // list of src.  The node is added at the begining
    int x=graph->array[dest].v_xcord;
    int y=graph->array[dest].v_ycord;
    struct Neigh_Node* newNode;
    if(type=='P'){
    	cnode++;
    	pnode++;
    	newNode = newNeigh_Node(dest,rname,dist,dist/Walk_Speed,dist/Car_Speed,dist/Public_Speed,t_den,x,y);
    }

    if(type=='C'){
    	cnode++;
    	 newNode = newNeigh_Node(dest,rname,dist,dist/Walk_Speed,dist/Car_Speed,INF,t_den,x,y);
    }

    if(type=='W'){
    	 newNode = newNeigh_Node(dest,rname,dist,dist/Walk_Speed,INF,INF,t_den,x,y);
    }

    newNode->next = graph->array[src].head;
    graph->array[src].head = newNode;

    // Since graph is undirected, add an edge from dest to src also
    x=graph->array[src].v_xcord;
    y=graph->array[src].v_ycord;
    if(type=='P')
        newNode = newNeigh_Node(src,rname,dist,dist/Walk_Speed,dist/Car_Speed,dist/Public_Speed,t_den,x,y);
    if(type=='C')
        newNode = newNeigh_Node(src,rname,dist,dist/Walk_Speed,dist/Car_Speed,INF,t_den,x,y);
    if(type=='W')
        newNode = newNeigh_Node(src,rname,dist,dist/Walk_Speed,INF,INF,t_den,x,y);
    newNode->next = graph->array[dest].head;
    graph->array[dest].head = newNode;
}


void printGraph(struct Graph* graph)
{
    int v;
    for (v = 0; v < graph->V; ++v)
    {
        struct Neigh_Node* temp = graph->array[v].head;
        cout<<"\n Adjacency list of vertex "<<graph->array[v].Place_Name<<" ("<<v<<")\n head ";
        while (temp)
	{
	    cout<<"-> "<<temp->dest;
	    cout<<" ("<<temp->Road_Name<<") ";
	    temp = temp->next;
	}
	cout<<"\n";
    }
}

struct dj_heap
{
    int  v;
    int dist;
};

struct heap
{
    int size;
    int cap;
    int *pos;
    struct dj_heap **array;
};


int checknodepresent_carpath(int n,struct Graph* graph);
int checknodepresent_publicpath(int n,struct Graph* graph);
void serach_adjcent_car(int n,struct Graph* graph);
void serach_typenode_adjcentnode(int s,struct Graph *graph,char type);
void serachadjcentnode(int,struct Graph* graph);
void printpath_dijstra(int,int);
void printpath_asearch(int arr[],int);
void printpath_asearch_public(int arr[],int);
void printpath_asearch(int arr[],int,int);
void path_viasomeplace_a_search(int s,int d,struct Graph* graph,int path[],int *pathsize,char type);
void path_viasomeplace_d_search(int s,int d,struct Graph* graph,int path[],int *pathsize,char type);
void path_car_graph_a_search(int ,int ,struct Graph*,int[],int*  );
void path_car_graph_d_search(int ,int ,struct Graph*,int[],int*  );
void path_public_graph_a_search(int ,int ,struct Graph*,int[],int* );
void path_nearestplace_a_search(int s,struct Graph* graph,int path[],int *pathsize,char type);
void path_nearestplace_d_search(int s,struct Graph* graph,int path[],int *pathsize,char type);
void path_public_graph_d_search(int s,int d,struct Graph* graph,int path[],int *pathsize);
void path_walk_graph_d_search(int s,int d,struct Graph* graph,int [],int*,char);
void path_walk_graph_a_search(int s,int d,struct Graph* graph,int [],int*,char);
void path_traffic_graph_a_search(int s,int d,struct Graph* graph,int path[],int *pathsize,char type);
void print_final_path(int path[],int pathsize);
void decrease_Key(struct heap* min_heap, int v, int dist);
bool isInheap(struct heap *min_heap, int v);
int isEmpty(struct heap* min_heap);
void min_heapify(struct heap* min_heap, int i);
void swap_djnode(struct dj_heap** a, struct dj_heap** b);
struct heap* createheap(int cap);
struct dj_heap* newdj_heap(int v,int dist);
void dijkstra_path(int,int,struct Graph *graph ,int dpath[],int *dsize,char);
void dijkstra(struct Graph* graph, int src,int prev[],char type);
void dijkstra_dist(struct Graph* graph, int src,char type,int dist[]);
void path_traffic_graph_d_search(int s,int d,struct Graph* graph,int path[],int *pathsize);




///////////////////to check weather node present in car path///////////////////////////////////

int checknodepresent_carpath(int n,struct Graph* graph)
{

		struct Neigh_Node* temp = graph->array[n].head;
        while (temp)
        {
            if(temp->c_time!=INF)
            return 1;
            else
            temp = temp->next;
        }
        return 0;

}

///////////////////to check weather node present in public path///////////////////////////////////

int checknodepresent_publicpath(int n,struct Graph* graph)
{

		struct Neigh_Node* temp = graph->array[n].head;
        while (temp)
        {
            if(temp->p_time!=INF)
            return 1;
            else
            temp = temp->next;
        }
        return 0;

}


///////////////////////////search for adjacent nodes present in car path/////////////////////////////////////////////
void serach_adjcent_car(int n,struct Graph* graph)
{
	struct Neigh_Node* temp = graph->array[n].head;
	int arr[10];
	int i=0;
	 while (temp)
        {
            if(temp->w_time!=INF)
			{
			arr[i]=temp->dest;
            i++;
            }
			temp = temp->next;
        }


for(int j=0;j<i;j++)
{
	if(visit[arr[j]]==-1)//node not visited
	{
	visit[arr[j]]=0;
	if(checknodepresent_carpath(arr[j],graph))
	{
		pathnodes[k]=arr[j];
		k++;


	}
	else
	{
	serach_adjcent_car(arr[j],graph);
	//	return ;
	}
}
}
}

void serachadjcentnode_car(int n,struct Graph* graph)//search for those adjacent nodes present in car path
{

	for(int i=0;i<V;i++)
	{
		visit[i]=-1;
		if(i==n)
		visit[i]=0;
	}
	k=0;
	serach_adjcent_car(n,graph);
}

///////////////////////////search for adjacent nodes present in public path/////////////////////////////////////////////
void serach_adjcent_public(int n,struct Graph* graph)
{
	struct Neigh_Node* temp = graph->array[n].head;
	int arr[10];
	int i=0;
	 while (temp)
	{
	    if(temp->w_time!=INF)
			{
			arr[i]=temp->dest;
	    i++;
	    }
			temp = temp->next;
	}


for(int j=0;j<i;j++)
{
	if(visit[arr[j]]==-1)//node not visited
	{
	visit[arr[j]]=0;
	if(checknodepresent_publicpath(arr[j],graph))
	{
		pathnodes[k]=arr[j];
		k++;


	}
	else
	{
	serach_adjcent_public(arr[j],graph);
	//	return ;
	}
}
}
}

void serachadjcentnode_public(int n,struct Graph* graph)//search for those adjacent nodes present in public path
{

	for(int i=0;i<V;i++)
	{
		visit[i]=-1;
		if(i==n)
		visit[i]=0;
	}
	k=0;
	serach_adjcent_public(n,graph);
}


/////////////////////////////////////////////////////////////////////
/////a* search algorthim ///////////////////////////////////////////



int car_weight(int k,int n,struct Graph *graph)
{
	struct Neigh_Node* temp = graph->array[n].head;
	while (temp)
	{
		if(temp->dest==k)
		return temp->c_time;
		temp=temp->next;
	}
}

int traffic_weight(int k,int n,struct Graph *graph)
{
	struct Neigh_Node* temp = graph->array[n].head;
	while (temp)
	{
		if(temp->dest==k)
		{//cout<<"traffic"<<temp->t_density;
			return temp->t_density;
		}

		temp=temp->next;
	}
}

int walk_weight(int k,int n,struct Graph *graph)
{
	struct Neigh_Node* temp = graph->array[n].head;
	while (temp)
	{
		if(temp->dest==k)
		return temp->w_time;
		temp=temp->next;
	}
}

int distance_weight(int k,int n,struct Graph *graph)
{
	struct Neigh_Node* temp = graph->array[n].head;
	while (temp)
	{
		if(temp->dest==k)
		return temp->distance;
		temp=temp->next;
	}
}

int public_weight(int k,int n,struct Graph *graph)
{
	struct Neigh_Node* temp = graph->array[n].head;
	while (temp)
	{
		if(temp->dest==k)
		return temp->p_time;
		temp=temp->next;
	}
}


//straight line distance
int h()
{
	return 0;
}

int h(int x1,int y1,int x2,int y2)
{
	return (sqrt((pow(x2-x1,2)+pow(y2-y1,2)))*50);
}

struct qnode{
	int fvalue;
	int gprev;
	int presentnode;
	int path[range];
	int psize;
}*data[range];

int heap_size=0;
int parent(int i) { return (i-1)/2; }
int left(int i) { return (2*i + 1); }
int right(int i) { return (2*i + 2); }


////////////////fuction for heap/////////////////////////////////
void swap(struct qnode *x, struct qnode *y)
{
    struct qnode temp = *x;
    *x = *y;
    *y = temp;
}



void insert_heap(struct qnode* q)
{

	if (heap_size == range)
    {
        cout << "\nOverflow\n";
        return;
    }

    heap_size++;
    int i = heap_size - 1;
   data[i] = q;

    while (i != 0 && data[parent(i)]->fvalue > data[i]->fvalue)
    {
       swap(data[i], data[parent(i)]);
       i = parent(i);
    }

}


void decreaseKey(int i, struct qnode* new_val)
{
    data[i] = new_val;
    while (i != 0 && data[parent(i)]->fvalue > data[i]->fvalue)
    {
       swap(data[i], data[parent(i)]);
       i = parent(i);
    }
}


void MinHeapify(int i)
{
    int l = left(i);
    int r = right(i);
    int smallest = i;
    if (l < heap_size && data[l]->fvalue < data[i]->fvalue)
        smallest = l;
    if (r < heap_size && data[r]->fvalue < data[smallest]->fvalue)
        smallest = r;
    if (smallest != i)
    {
        swap(data[i], data[smallest]);
        MinHeapify(smallest);
    }
}
struct qnode* extractMin()
{
    if (heap_size <= 0)
        return NULL;
    if (heap_size == 1)
    {
        heap_size--;
        return data[0];
    }

    struct qnode * root = data[0];
    data[0] = data[heap_size-1];
    heap_size--;
    MinHeapify(0);

    return root;
}


void deleteKey(int i)
{
	struct qnode *temp=(struct qnode*) malloc(sizeof(struct qnode));
	temp->fvalue=range;
	temp->gprev=range;
	temp->presentnode=range;
	temp->path[0]=range;
	temp->psize=0;;
    decreaseKey(i,temp);
    extractMin();
}

///////////////////heap function end /////////////////////////////////////



void insertqueue(struct Graph *graph,int pn,int k,int gprev,int arr[],int psize,int fvalue,char type)
{
	//cout<<"\n inside insert";

		struct qnode *temp=(struct qnode*) malloc(sizeof(struct qnode));
		temp->presentnode=pn;
		//cout<<"\n"<<adjnodes[j];
		temp->fvalue=fvalue;
		if(type=='c')
		{
			temp->gprev=car_weight(pn,k,graph)+gprev;
		}
		if(type=='w')
		{
			temp->gprev=walk_weight(pn,k,graph)+gprev;
		}
		if(type=='p')
		{
			temp->gprev=public_weight(pn,k,graph)+gprev;
		}
		if(type=='d')
		{
			temp->gprev=distance_weight(pn,k,graph)+gprev;
		}
		if(type=='t')
		{
			temp->gprev=traffic_weight(pn,k,graph)+gprev;
		}
		int i;

				for( i=0;i<psize;i++)
				{
					temp->path[i]=arr[i];

				}

		temp->path[i]=pn;
		temp->psize=i+1;
		insert_heap(temp);
}


void serach_adjcent(int n,struct Graph* graph,char type)
{
	//cout<<"\n search adjaenct node";
	struct Neigh_Node* temp = graph->array[n].head;
	//int arr[10];
	int i=0;
	if(type=='w')
	{
 	 while (temp)
        {
			adjnodes[i]=temp->dest;
            i++;

			temp = temp->next;
        }

 			 sizeadjnodes=i;
	}
	if(type=='p')
	{
 	while (temp)
        {
        	if(temp->p_time!=INF)
        	{

			adjnodes[i]=temp->dest;
            i++;
       		 }
			temp = temp->next;
        }

	sizeadjnodes=i;
	}

 	if(type='c')
 	{
 	 while (temp)
        {
        	if(temp->c_time!=INF)
			{
			adjnodes[i]=temp->dest;
            i++;
        	}
			temp = temp->next;
        }

	sizeadjnodes=i;
 	}
 	 	if(type='t')
 	{
 	 while (temp)
        {
        	if(temp->c_time!=INF)
			{
			adjnodes[i]=temp->dest;
            i++;
        	}
			temp = temp->next;
        }

	sizeadjnodes=i;
 	}

}

void serach_typenode_adjcentnode(int n,struct Graph* graph,char type,int nodes[],int *size)
{
	//cout<<"\n search adjaenct node";
	char t ;

	int i=0,j=1;
		while (j<vertex)
        {
        	t = graph->array[j].Place_Type;
        	//cout<<"\t"<<t;
			if(t==type)
        	{
			nodes[i]=j;
			i++;
       		 }
			j++;
        }
	*size=i;
}


struct qnode* a_search(int src,int des,struct Graph *graph,char type)
{

	struct qnode * start=(struct qnode*) malloc(sizeof(struct qnode));
	struct qnode *p;

	heap_size=0;

	start->fvalue=0+h();
	start->presentnode=src;
	start->gprev=0;
	start->path[0]=src;
	start->psize=1;
	insert_heap(start);

	int visit[vertex];

	for(int j=0;j<vertex;j++)
	{
		visit[j]=-1;
	}

	int n;

	do
	{
		p=extractMin();
		visit[p->presentnode]=0;
		n=p->presentnode;//p->present node not working directly

		if(n==des)
		return p;

		if(type=='c'||type=='C'||type=='t'||type=='T')
		{
				serach_adjcent(p->presentnode,graph,'c');

				for(int j=0;j<sizeadjnodes;j++)
				{

				int i=0;
				int x1,x2,y1,y2;
				x1=graph->array[adjnodes[j]].v_xcord;
				y1=graph->array[adjnodes[j]].v_ycord;
				x2=graph->array[des].v_xcord;
				y2=graph->array[des].v_ycord;
				int fvalue;
				if(type=='c')
					fvalue=car_weight(adjnodes[j],p->presentnode,graph)+p->gprev+h(x1,y1,x2,y2);
					if(type=='t')
					fvalue=traffic_weight(adjnodes[j],p->presentnode,graph)+p->gprev+h(x1,y1,x2,y2);

				//cout<<"\n node deleted "<<p->presentnode<<"  "<<p->fvalue<<" "<<p->gprev;
				if(visit[adjnodes[j]]==-1 || fvalue<p->fvalue)
				insertqueue(graph,adjnodes[j],p->presentnode,p->gprev,p->path,p->psize,fvalue,type);
				}
		}

		if(type=='w'||type=='W')
		{

				serach_adjcent(p->presentnode,graph,'w');

				for(int j=0;j<sizeadjnodes;j++)
				{

				int i=0;
				int x1,x2,y1,y2;
				x1=graph->array[adjnodes[j]].v_xcord;
				y1=graph->array[adjnodes[j]].v_ycord;
				x2=graph->array[des].v_xcord;
				y2=graph->array[des].v_ycord;
				int fvalue;

					fvalue=walk_weight(adjnodes[j],p->presentnode,graph)+p->gprev+h(x1,y1,x2,y2);

				//cout<<"\n node deleted "<<p->presentnode<<"  "<<p->fvalue<<" "<<p->gprev;
				if(visit[adjnodes[j]]==-1 || fvalue<p->fvalue)
				insertqueue(graph,adjnodes[j],p->presentnode,p->gprev,p->path,p->psize,fvalue,type);
				}
		}

		if(type=='p'||type=='P')
		{

				serach_adjcent(p->presentnode,graph,'p');

				for(int j=0;j<sizeadjnodes;j++)
				{

				int i=0;
				int x1,x2,y1,y2;
				x1=graph->array[adjnodes[j]].v_xcord;
				y1=graph->array[adjnodes[j]].v_ycord;
				x2=graph->array[des].v_xcord;
				y2=graph->array[des].v_ycord;
				int fvalue;

					fvalue=public_weight(adjnodes[j],p->presentnode,graph)+p->gprev+h(x1,y1,x2,y2);

				//cout<<"\n node deleted "<<p->presentnode<<"  "<<p->fvalue<<" "<<p->gprev;
				if(visit[adjnodes[j]]==-1 || fvalue<p->fvalue)
				insertqueue(graph,adjnodes[j],p->presentnode,p->gprev,p->path,p->psize,fvalue,type);
				}

		}
		if(type=='d'||type=='D')
		{

				serach_adjcent(p->presentnode,graph,'w');

				for(int j=0;j<sizeadjnodes;j++)
				{

				int i=0;
				int x1,x2,y1,y2;
				x1=graph->array[adjnodes[j]].v_xcord;
				y1=graph->array[adjnodes[j]].v_ycord;
				x2=graph->array[des].v_xcord;
				y2=graph->array[des].v_ycord;
				int fvalue;

					fvalue=distance_weight(adjnodes[j],p->presentnode,graph)+p->gprev+h(x1,y1,x2,y2);

				//cout<<"\n node deleted "<<p->presentnode<<"  "<<p->fvalue<<" "<<p->gprev;
				if(visit[adjnodes[j]]==-1 || fvalue<p->fvalue)
				insertqueue(graph,adjnodes[j],p->presentnode,p->gprev,p->path,p->psize,fvalue,type);
				}

		}
	}while(1);
	return p;
}


//////////////////////////////graphics///////////////////////////////////
/*
void init_graphics()      //initialize graphics driver
{
   // request auto detection
   int gdriver = DETECT, gmode, errorcode;

   // initialize graphics and local variables
   initgraph(&gdriver, &gmode, "C:\\TURBOC3\\BGI");

   // read result of initialization
   errorcode = graphresult();
   if (errorcode != grOk)
   // an error occurred
   {
     printf("Graphics error: %s\n", grapherrormsg(errorcode));
     printf("Press any key to halt:");
     getch();
     // terminate with an error code
     exit(1);
   }
}  */

void Draw_Graph(Graph* graph)
{
    cleardevice();
    setcolor(YELLOW);
    settextstyle(2,0,9);
    settextjustify(1,2);
    outtextxy(320*SCALE,15*SCALE,"LOS SANTOS MAP");
    settextjustify(0,2);
    for (int v = 1; v <= graph->V; ++v)
    {
        int vx=graph->array[v].v_xcord;
        int vy=graph->array[v].v_ycord;
        char vtype=graph->array[v].Place_Type;
        char* string=graph->array[v].Place_Name;
        itoa(v,string,10);

        Place_Vertex(vx,vy,vtype);

        struct Neigh_Node* temp = graph->array[v].head;
        int nx,ny;
        char ntype;
        int color;
        while (temp)
        {
            nx=temp->n_xcord;
            ny=temp->n_ycord;

            if(temp->p_time!=INF)
                {
                    ntype='P';
                    color=WHITE;
                }
            else if(temp->c_time!=INF)
                {
                    ntype='C';
                    color=WHITE;
                }

            else
                {
                    ntype='W';
                    color=WHITE;
                }

            Draw_Road(vx,vy,nx,ny,ntype,color);
            Display_Name(vx,vy+2,string,BLUE);

	    temp = temp->next;
	}
    }
    setcolor(WHITE);
    setlinestyle(3,1,1);
    line(500*SCALE,100*SCALE,525*SCALE,100*SCALE);
    setcolor(WHITE);
    settextstyle(2,0,5);
    outtextxy(535*SCALE,100*SCALE,"Walking Path");

    setcolor(WHITE);
    setlinestyle(0,1,1);
    line(500*SCALE,120*SCALE,525*SCALE,120*SCALE);
    setcolor(WHITE);
    settextstyle(2,0,5);
    outtextxy(535*SCALE,120*SCALE,"Car Route");

    setcolor(WHITE);
    setlinestyle(0,1,3);
    line(500*SCALE,140*SCALE,525*SCALE,140*SCALE);
    setcolor(WHITE);
    settextstyle(2,0,5);
    outtextxy(535*SCALE,140*SCALE,"Public Transport Route");


}

void Place_Vertex(int x, int y, char type)
{
    circle(x*SCALE,y*SCALE,2);
}

void Draw_Road(int x1, int y1, int x2, int y2, char type,int color)
{
    if(type=='W')
    {
        setcolor(color);
        setlinestyle(1,1,1);
        line(x1*SCALE,y1*SCALE,x2*SCALE,y2*SCALE);
    }
    else if(type=='C')
    {
        setcolor(color);
        setlinestyle(0,1,1);
        line(x1*SCALE,y1*SCALE,x2*SCALE,y2*SCALE);
    }
    else if(type=='P')
    {
        setcolor(color);
        setlinestyle(0,1,3);
        line(x1*SCALE,y1*SCALE,x2*SCALE,y2*SCALE);
    }

    else
    {
        setcolor(color);
        setlinestyle(0,1,1);
        line(x1*SCALE,y1*SCALE,x2*SCALE,y2*SCALE);
    }


}

void Display_Name(int x, int y, char* string,int color)
{
    setcolor(color);
    settextstyle(2,0,4);
    outtextxy(x*SCALE,y*SCALE,string);
}


void Display_Path(Graph* graph,int path[],int path_len,int color)
{
    for(int i=0;i<path_len-1;i++)
    {
	int x1=graph->array[path[i]].v_xcord;
	int y1=graph->array[path[i]].v_ycord;
	char type1=graph->array[path[i]].Place_Type;
	int x2=graph->array[path[i+1]].v_xcord;
	int y2=graph->array[path[i+1]].v_ycord;
	int type2=graph->array[path[i+1]].Place_Type;
	char* string=graph->array[i].Place_Name;
	char ntype;

	//setcolor(1);
	Place_Vertex(x1,y1,type1);
	delay(100);

	    struct Neigh_Node* temp=graph->array[path[i]].head;
        while(temp->dest!=path[i+1])
            temp=temp->next;
        if(temp->p_time!=INF)
                    ntype='P';

        else if(temp->c_time!=INF)
                    ntype='C';

        else
                    ntype='W';

	Draw_Road(x1,y1,x2,y2,ntype,color);
	delay(100);
	Place_Vertex(x2,y2,type2);
	itoa(path[i],string,10);
	Display_Name(x1,y1+2,string,BLUE);
	itoa(path[i+1],string,10);
    Display_Name(x2,y2+2,string,BLUE);

    }

}

void dijkstra_graph(Graph* graph,int path[],int path_len)
{

    //Draw_Graph(graph);
    setcolor(RED);
    setlinestyle(0,1,1);
    line(200*SCALE,430*SCALE,225*SCALE,430*SCALE);
    setcolor(WHITE);
    settextstyle(2,0,5);
    outtextxy(235*SCALE,430*SCALE,"dijkstra Algorithm Path");
    Display_Path(graph,path,path_len,RED);
}

void astar_graph(Graph* graph,int path[],int path_len)
{

    Draw_Graph(graph);
    setcolor(GREEN);
    setlinestyle(0,1,1);
    line(325*SCALE,430*SCALE,350*SCALE,430*SCALE);
    setcolor(WHITE);
    settextstyle(2,0,5);
    outtextxy(360*SCALE,430*SCALE,"A* Algorithm Path");
    Display_Path(graph,path,path_len,GREEN);
}


void welcome()
{
    int i,j;
    setcolor(WHITE);
    settextjustify(1,0);
    settextstyle(3,0,8);
    outtextxy(320*SCALE,150*SCALE,"FAST");
    outtextxy(320*SCALE,230*SCALE,"AND");
    outtextxy(320*SCALE,310*SCALE,"FURIOUS");
    settextstyle(2,0,6);
    setcolor(BLUE);
    outtextxy(325*SCALE,370*SCALE,"Loading");
    for(j=1;j<=3;j++)
    {
        for(i=0;i<=3;i++)
        {
            putpixel((360+(i*10))*SCALE,370*SCALE,WHITE);
            putpixel((288-(i*10))*SCALE,370*SCALE,WHITE);
            delay(250);
        }
        for(i=0;i<=3;i++)
        {
            putpixel((390-(i*10))*SCALE,370*SCALE,BLACK);
            putpixel((258+(i*10))*SCALE,370*SCALE,BLACK);
            delay(250);
        }
    }
    cleardevice();
    settextjustify(0,2);     //Default Text Justification
}

////////////////////////////////////////////////////////////////////////////



int main()
{
    // create the graph given
    int s,d,flag1,flag2;
    int V =111;
    struct Graph* graph = createGraph(V);
    Add_Roads(graph);

   // printGraph(graph);
   // getch();
   // clrscr();
  //  init_graphics();
    initwindow(960,720);
    welcome();
    Draw_Graph(graph);
    delay(1000);
    //int path[]={47,48,49,50,46,4,5,6,9,11,12,13,14,15,108,43,41,42,67,68};
   // getch();*/
  //  return 0;
     //clrscr();
	// print the adjacency list representation of the above graph
    //printGraph(graph);
    char ch='y';
    int choice1,choice2,choice3;
	int path[MAX];
	int pathsize;


    //clrscr();
	while(ch=='y'||ch=='Y')
    {


	cout<<"\n1.Simple path";
	cout<<"\n2.Path Via some other palce";
	cout<<"\n3.Search for nearest hotel /pterol pump etc";
	cout<<"\nenter your choice";
	cin>>choice1;
	switch(choice1)
		{
			case 1:	do
			{
					cout<<"\nSimple path";
					cout<<"\n-----------";
					cout<<"\nenter source destination";
					cin>>s>>d;
					cout<<"\nGet path according to:";
					cout<<"\n1. Traffic";
					cout<<"\n2. Distance";
					cout<<"\n3. Time";
					cout<<"\nEnter your choice";
					cin>>choice2;
					switch(choice2)
					{
						case 1:
								path_traffic_graph_a_search(s,d,graph,path,&pathsize,'t');
								cout<<"\nA* search path:";
								cout<<"\n----------------";
								print_final_path(path,pathsize);
								astar_graph(graph,path,pathsize);
								//getch() ;

							//	path_traffic_graph_d_search(s,d,graph,path,&pathsize);
							 //	dijkstra_path(s,d,graph,path,&pathsize,'t');
							//	cout<<"\nDijkstra path:";
							//	cout<<"\n----------------";
							//	print_final_path(path,pathsize);

							break;
						case 2:path_walk_graph_a_search(s,d,graph,path,&pathsize,'d');
							cout<<"\n\nA* search path:";
							print_final_path(path,pathsize);
							astar_graph(graph,path,pathsize);

						 	dijkstra_path(s,d,graph,path,&pathsize,'d');
							cout<<"\n\nDijkstra path:";
							cout<<"\n---------------";
							print_final_path(path,pathsize);
							dijkstra_graph(graph,path,pathsize);

							break;
						case 3:
							do
							{

							cout<<"\nGET DIRECTION";
						cout<<"\n-------------";
						cout<<"\n1.Walking path";
						cout<<"\n2.Car path";
						cout<<"\n3.Public path";
						cout<<"\nEnter your choice";
							cin>>choice3;

							switch(choice3)
								{
									case 1:
											path_walk_graph_a_search(s,d,graph,path,&pathsize,'w');
											cout<<"\nA* search path:";
											cout<<"\n---------------";
											print_final_path(path,pathsize);
											astar_graph(graph,path,pathsize);
											getch();

										 	dijkstra_path(s,d,graph,path,&pathsize,'w');
											cout<<"\n\nDijkstra path:";
											cout<<"\n---------------";
											print_final_path(path,pathsize);
											dijkstra_graph(graph,path,pathsize);

										break;
									case 2:
											path_car_graph_a_search(s,d,graph,path,&pathsize);
											cout<<"\n\nA* search path:";
											cout<<"\n---------------";
											print_final_path(path,pathsize);
											astar_graph(graph,path,pathsize);
											getch() ;

										 	//path_car_graph_d_search(s,d,graph,path,&pathsize);
											cout<<"\n\nDijkstra path:";
											cout<<"\n---------------";
											print_final_path(path,pathsize);
											dijkstra_graph(graph,path,pathsize);

										break;
									case 3:
											path_public_graph_a_search(s,d,graph,path,&pathsize);
											cout<<"\n\nA* search path:";
											cout<<"\n---------------";
											print_final_path(path,pathsize);
											astar_graph(graph,path,pathsize);
											getch() ;

										 	path_public_graph_d_search(s,d,graph,path,&pathsize);
											cout<<"\n\nDijkstra path:";
											cout<<"\n---------------";
											print_final_path(path,pathsize);
											dijkstra_graph(graph,path,pathsize);

										break;
									default:cout<<"\nwrong choice";
										break;
								}
							cout<<"\nDo you want to continue......";
							cin>>ch;
							}while(ch=='y'||ch=='Y');
							break;
						default:cout<<"\nwrong choice";
								break;
					}
			cout<<"\nDo you want to continue......";
			cin>>ch;

			}while(ch=='y'||ch=='Y');
					break;
		    case 2:	do
			{
					cout<<"\nPath Via some other palce";
					cout<<"\n--------------------------";
					cout<<"\nenter source destination";
					cin>>s>>d;
					cout<<"\nWant to go via \n1.hotel\n2.petrol pump\n3.shopping mall";
					cout<<"\nEnter your choice";
					cin>>choice2;
					switch(choice2)
					{
						case 1:path_viasomeplace_a_search(s,d,graph,path,&pathsize,'H');
							cout<<"\n\nA* search path:";
							cout<<"\n---------------";
							print_final_path(path,pathsize);
							astar_graph(graph,path,pathsize);
							getch() ;

						 	path_viasomeplace_d_search(s,d,graph,path,&pathsize,'H');
							cout<<"\n\nDijkstra path:";
							cout<<"\n---------------";
							print_final_path(path,pathsize);
							dijkstra_graph(graph,path,pathsize);

							break;
						case 2:path_viasomeplace_a_search(s,d,graph,path,&pathsize,'P');
							cout<<"\nA* search path:";
							cout<<"\n---------------";
							print_final_path(path,pathsize);
							astar_graph(graph,path,pathsize);
							getch() ;

						 	path_viasomeplace_d_search(s,d,graph,path,&pathsize,'P');
							cout<<"\n\nDijkstra path:";
							cout<<"\n---------------";
							print_final_path(path,pathsize);
							dijkstra_graph(graph,path,pathsize);

							break;
						case 3:path_viasomeplace_a_search(s,d,graph,path,&pathsize,'S');
							cout<<"\nA* search path:";
							cout<<"\n---------------";
							print_final_path(path,pathsize);
							astar_graph(graph,path,pathsize);
							getch() ;

						 	path_viasomeplace_d_search(s,d,graph,path,&pathsize,'S');
							cout<<"\n\nDijkstra path:";
							cout<<"\n---------------";
							print_final_path(path,pathsize);
							dijkstra_graph(graph,path,pathsize);

							break;
						default:cout<<"\nwrong choice";
								break;
					}
			cout<<"\nDo you want to continue......";
			cin>>ch;

			}while(ch=='y'||ch=='Y');
					break;
			case 3:	do
			{
					cout<<"\nSearch for nearest hotel /pterol pump etc";
					cout<<"\n------------------------------------------";
					cout<<"\nEnter source ";
					cin>>s;
					cout<<"\nWant to search for : \n1.hotel\n2.petrol pump\n3.shopping mall";
					cout<<"Enter your choice";
					cin>>choice2;
					switch(choice2)
					{
						case 1:path_nearestplace_a_search(s,graph,path,&pathsize,'H');
							cout<<"\nA* search path:";
							cout<<"\n---------------";
							print_final_path(path,pathsize);
							astar_graph(graph,path,pathsize);
							getch() ;

						 	path_nearestplace_d_search(s,graph,path,&pathsize,'H');
							cout<<"\n\nDijkstra path:";
							cout<<"\n---------------";
							print_final_path(path,pathsize);
							dijkstra_graph(graph,path,pathsize);

							break;
						case 2:path_nearestplace_a_search(s,graph,path,&pathsize,'P');
							cout<<"\nA* search path:";
							cout<<"\n---------------";
							print_final_path(path,pathsize);
							astar_graph(graph,path,pathsize);
							getch() ;

						 	path_nearestplace_d_search(s,graph,path,&pathsize,'P');
							cout<<"\n\nDijkstra path:";
							cout<<"\n---------------";
							print_final_path(path,pathsize);
							dijkstra_graph(graph,path,pathsize);

							break;
						case 3:
							path_nearestplace_a_search(s,graph,path,&pathsize,'S');
							cout<<"\nA* search path:";
							cout<<"\n---------------";
							print_final_path(path,pathsize);
							astar_graph(graph,path,pathsize);
							getch() ;

						 	path_nearestplace_d_search(s,graph,path,&pathsize,'S');
							cout<<"\n\nDijkstra path:";
							cout<<"\n---------------";
							print_final_path(path,pathsize);
							dijkstra_graph(graph,path,pathsize);

							break;
						default:cout<<"\nwrong choice";
								break;
					}
					cout<<"\nDo you want to continue......";
			cin>>ch;

			}while(ch=='y'||ch=='Y');
			break;
			default:cout<<"\nwrong choice";
								break;
		}
		cout<<"\nDo you want to continue......";
		cin>>ch;
	}

	//clrscr();
	//getch();
	return 0;
}

struct qnode * intalize_min()
{
struct qnode *min=(struct qnode*) malloc(sizeof(struct qnode));
	min->gprev=range;
 	min->fvalue=range;
 	min->presentnode=range;
 	min->path[0]=range;
 	min->psize=range;
	 return min;
}

void path_nearestplace_a_search(int s,struct Graph* graph,int path[],int *pathsize,char type)
{//cout<<"hi";
		char *t;
		if(type=='H'){t = new char[strlen("HOTEL")+1];
		strcpy(t,"HOTEL");
		}
		if(type=='P'){t = new char[strlen("PETROL PUMP")+1];
		strcpy(t,"PETROL PUMP");
		}
		int node[MAX];
		int size;
		struct qnode *min,*p;/*
		if(graph->array[s].Place_Type==type)
		{
			cout<<"\n\n the nearest "<<t<<" is "<<s;
			path[0]=s;
			*pathsize=1;
			return ;
		}*/
		min=intalize_min();
		int i;

		serach_typenode_adjcentnode(s,graph,type,node,&size);
		for(i=0;i<size;i++)
			 {
			 //	cout<<"\nadjnode: "<<node[i];
	 			p=a_search(s,node[i],graph,'w');
			 	if(min->gprev>p->gprev)
			 		{
			 			min=p;
			 		}
			 }
			 cout<<"\n\n the nearest "<<t<<" is "<<min->path[min->psize-1];
		//	 printpath_asearch(min->path,min->psize,0);
		for(i=0;i<min->psize;i++)
		{
			path[i]=min->path[i];
		}
		*pathsize=i;
}

void path_nearestplace_d_search(int s,struct Graph* graph,int path[],int *pathsize,char type)
{
        	char *t;
		if(type=='H'){t = new char[strlen("HOTEL")+1];
		strcpy(t,"HOTEL");
		}
		if(type=='P'){t = new char[strlen("PETROL PUMP")+1];
		strcpy(t,"PETROL PUMP");
		}

		/*if(graph->array[s].Place_Type==type)
		{
		cout<<"\n\n the nearest "<<t<<" is "<<s;
			path[0]=s;
			*pathsize=1;
			return;
		}*/
		int node[MAX];
		int size;
		int i,cnode,min=INT_MAX,dp[MAX],dsize;
		int dist[vertex];
		dijkstra_dist(graph,s,'w',dist);
	//	cout<<"\nhi";
		serach_typenode_adjcentnode(s,graph,type,node,&size);
		for(i=0;i<size;i++)
		 {
		// cout<<"\nnode "<<node[i]<<" dist "<<dist[node[i]];
			 if(dist[node[i]]<=min)
		 	{
		 		min=dist[node[i]];
		 		cnode=node[i];
		 	}
		 //	cout<<"cnode"<<cnode;
		 }
		 //cout<<"cnode"<<cnode;
		 cout<<"\n\n the nearest "<<t<<" is "<<cnode;
	 	 dijkstra_path(s,cnode,graph,dp,&dsize,'w');
	 	 //cout<<"hi 2";
		for(i=0;i<dsize;i++)
		{
			path[i]=dp[i];
		}
		*pathsize=i;
}

void path_viasomeplace_d_search(int s,int d,struct Graph* graph,int path[],int *pathsize,char type)
{
 	 		 	char *t;
		if(type=='H'){t = new char[strlen("HOTEL")+1];
		strcpy(t,"HOTEL");
		}
		if(type=='P'){t = new char[strlen("PETROL PUMP")+1];
		strcpy(t,"PETROL PUMP");
		}
		int i,size,node[MAX];
		int cnode,min=INT_MAX,dp[MAX],dsize;
		int dist[vertex];
		if(graph->array[d].Place_Type==type)
		{
			cout<<"\n\n the nearest "<<t<<" is "<<d;
			dijkstra_path(s,d,graph,dp,&dsize,'w');
			int j=0;
			for(i=0;i<dsize;i++)
			{
				path[j]=dp[i];
				j++;
			}
				*pathsize=j;
			return;
		}

		dijkstra_dist(graph,s,'w',dist);
		serach_typenode_adjcentnode(s,graph,type,node,&size);
		for(i=0;i<size;i++)
			 {
	 		 if(dist[node[i]]<min)
			 	{
			 		min=dist[node[i]];
			 		cnode=node[i];
			 	}
			 }

		cout<<"\n\n the nearest "<<t<<" is "<<cnode;
		dijkstra_path(s,cnode,graph,dp,&dsize,'w');
		int j=0;
		for(i=0;i<dsize;i++)
		{
			path[j]=dp[i];
			j++;
		}
		dijkstra_path(cnode,d,graph,dp,&dsize,'w');
		for(i=1;i<dsize;i++)
		{
			path[j]=dp[i];
			j++;
		}
		*pathsize=j;
}



void path_viasomeplace_a_search(int s,int d,struct Graph* graph,int path[],int *pathsize,char type)
{
		struct qnode *min,*p;
		min=intalize_min();
		int i,j,size,node[MAX];
			char *t;
		if(type=='H'){t = new char[strlen("HOTEL")+1];
		strcpy(t,"HOTEL");
		}
		if(type=='P'){t = new char[strlen("PETROL PUMP")+1];
		strcpy(t,"PETROL PUMP");
		}

		/*if(graph->array[d].Place_Type==type)
		{
			cout<<"\ndestination is palce of that type u r searching for";
			p=a_search(s,d,graph,'w');
			for(i=0;i<p->psize;i++)
			{
				path[j]=p->path[i];
				j++;
			}
			*pathsize=j;
			return;
		}*/

		serach_typenode_adjcentnode(s,graph,type,node,&size);
		for(i=0;i<size;i++)
			 {
	 			p=a_search(s,node[i],graph,'w');
			 	if(min->gprev>p->gprev)
			 		{
			 			min=p;
			 		}
			 }

		p=a_search(min->path[min->psize-1],d,graph,'w');
		cout<<"\n\n the nearest "<<t<<" is "<<min->path[min->psize-1]<<"\n";
		//printpath_asearch(min->path,min->psize,0);
 		//printpath_asearch(p->path,p->psize,0);
		j=0;
		for(i=0;i<min->psize;i++)
		{
			path[j]=min->path[i];
			cout<<"\t"<<min->path[i];
			j++;
		}
		for(i=1;i<p->psize;i++)
		{
			path[j]=p->path[i];
			j++;
		}
		*pathsize=j;
}

void path_traffic_graph_a_search(int s,int d,struct Graph* graph,int path[],int *pathsize,char type)
{

    int flag1=checknodepresent_carpath(s,graph);
	int flag2=checknodepresent_carpath(d,graph);



	if(flag1==1&&flag2==1)
		{
		struct qnode *min;
			min=a_search(s,d,graph,type);
			int i;
			for(i=0;i<min->psize;i++)
			{
				path[i]=min->path[i];
			}
			*pathsize=i;
			}
			else
			{
			cout<<"\nnodes are not on traffic path.";
			*pathsize=-1;
			}
 		//printpath_asearch(min->path,min->psize,0);
}
void path_walk_graph_a_search(int s,int d,struct Graph* graph,int path[],int *pathsize,char type)
{

		struct qnode *min;
			min=a_search(s,d,graph,type);
			int i;
			for(i=0;i<min->psize;i++)
			{
				path[i]=min->path[i];
			}
			*pathsize=i;

 		//printpath_asearch(min->path,min->psize,0);
}

int path_walk_car_walk(int s,int d,struct Graph* graph,int path[])
{
		k=0;
		int i;
		struct qnode *walk_car,*car_car,*car_walk;
		struct qnode *min=intalize_min();

	 	serachadjcentnode_car(s,graph);//carpathnodes[]
		for(i=0;i<k;i++)
			 {
	 			walk_car=a_search(s,pathnodes[i],graph,'w');
			 	if(min->gprev>walk_car->gprev)
			 		{
			 			min=walk_car;
			 		}
			 }
			 walk_car=min;

		min=intalize_min();

		serachadjcentnode_car(d,graph);//carpathnodes[]
		for(i=0;i<k;i++)
			 {
	 			car_walk=a_search(pathnodes[i],d,graph,'w');
			 	if(min->gprev>car_walk->gprev)
			 		{
			 			min=car_walk;
			 		}
			 }
			 car_walk=min;

		car_car=a_search(walk_car->path[walk_car->psize-1],car_walk->path[0],graph,'c');
//		printpath_asearch(walk_car->path,walk_car->psize,0);
//		cout<<"\nthe path by car:";
//		printpath_asearch(car_car->path,car_car->psize);
//		printpath_asearch(car_walk->path,car_walk->psize,0);
		int j;
		for( j=0;j<walk_car->psize;j++)
		{
			path[j]=walk_car->path[j];
		}

		for(i=1;i<car_car->psize;i++)
		{
			path[j]=car_car->path[i];
			j++;
		}

		for(i=1;i<car_walk->psize;i++)
		{
			path[j]=car_walk->path[i];
			j++;
		}
		return j;

}

void path_car_graph_a_search(int s,int d,struct Graph* graph,int path[],int *pathsize)
{
	struct qnode *p;
	int i;
	struct qnode *min=(struct qnode*) malloc(sizeof(struct qnode));
	int flag1=checknodepresent_carpath(s,graph);
	int flag2=checknodepresent_carpath(d,graph);

 	min->gprev=range;
 	min->fvalue=range;
 	min->presentnode=range;
 	min->path[0]=range;
 	min->psize=range;

 	if(flag1==0&&flag2==0)
 	{
 		*pathsize=path_walk_car_walk(s,d,graph,path);
 	}
 	if(flag1==1&&flag2==1)
 	{
 		min=a_search(s,d,graph,'c');
		for(i=0;i<min->psize;i++)
		{
			path[i]=min->path[i];
		}
		*pathsize=i;
 		//cout<<"\nthe path by car:";
 		//printpath_asearch(min->path,min->psize);
 	}

 	if((flag1==1&&flag2==0)||(flag1==0&&flag2==1))
 	{
 		//find adjacent node to one not in graph which exsist in graph
	 	if(flag1==0)
	 	{
	 		k=0;
	 		serachadjcentnode_car(s,graph);//carpathnodes[]
	 		for(i=0;i<k;i++)
			 {
			 	p=a_search(s,pathnodes[i],graph,'w');

			 	if(min->gprev>p->gprev)
			 		{
			 			min=p;
			 		}
			 }
		//	printpath_asearch(min->path,min->psize,0);
			p=a_search(min->path[min->psize-1],d,graph,'c');
		//	cout<<"\nthe path by car:";
		//	printpath_asearch(p->path,p->psize);
			int j=0;
			for(i=0;i<min->psize;i++)
				{
					path[j]=min->path[i];
					j++;
				}
			for(i=1;i<p->psize;i++)
			{
				path[j]=p->path[i];
				j++;
			}
			*pathsize=j;
	 	}
	 	if(flag2==0)
	 	{	k=0;
	 		serachadjcentnode_car(d,graph);
	 		for(i=0;i<k;i++)
			 {
	 			p=a_search(pathnodes[i],d,graph,'w');
			 	if(min->gprev>p->gprev)
			 		{
			 			min=p;
			 		}
			 }
			p=a_search(s,min->path[0],graph,'c');
			//cout<<"\nthe path by car:";
			//printpath_asearch(p->path,p->psize);
			//printpath_asearch(min->path,min->psize,0);
			int j=0;
			for(i=0;i<p->psize;i++)
			{
				path[j]=p->path[i];
				j++;
			}
			for(i=1;i<min->psize;i++)
				{
					path[j]=min->path[i];
					j++;
				}

			*pathsize=j;
	 	}
 	}
}


int path_walk_car_walk_d(int s,int d,struct Graph* graph,int path[])
{
		k=0;
		int i;
 		int min=INT_MAX;
 		int carnode,carnode2,dp[MAX],dsize,dp2[MAX],dsize2;
 		serachadjcentnode_car(s,graph);//carpathnodes[]
 			int dist[vertex];
	 	dijkstra_dist(graph,s,'w',dist);
		 for(i=0;i<k;i++)
		 {
			 if(dist[pathnodes[i]]<min)
		 	{
		 		min=dist[pathnodes[i]];
		 		carnode=pathnodes[i];
		 	}
		 }
	 	 dijkstra_path(s,carnode,graph,dp,&dsize,'w');
 		int j;
		for( j=0;j<dsize;j++)
		{
			path[j]=dp[j];
		}
		min=INT_MAX;
		serachadjcentnode_car(d,graph);//carpathnodes[]
	 	dijkstra_dist(graph,d,'w',dist);
		 for(i=0;i<k;i++)
		 {
			 if(dist[pathnodes[i]]<min)
		 	{
		 		min=dist[pathnodes[i]];
		 		carnode2=pathnodes[i];
		 	}
		 }
	 	 dijkstra_path(carnode2,d,graph,dp2,&dsize2,'w');
	 	 dijkstra_path(carnode,carnode2,graph,dp,&dsize,'c');
		for(i=1;i<dsize;i++)
		{
			path[j]=dp[i];
			j++;
		}

		for(i=1;i<dsize2;i++)
		{
			path[j]=dp2[i];
			j++;
		}
		return j;

}


void path_car_graph_d_search(int s,int d,struct Graph* graph,int path[],int *pathsize)
{
	int i;
	int flag1=checknodepresent_carpath(s,graph);
	int flag2=checknodepresent_carpath(d,graph);
	int dp[MAX],dsize;
		int dist[vertex];
 	if(flag1==0&&flag2==0)
 	{
        path_walk_car_walk_d(s,d,graph,path);
 	}
 	if(flag1==1&&flag2==1)
 	{
 		dijkstra_path(s,d,graph,path,&dsize,'c');
		 *pathsize=dsize;
 	}

 	if((flag1==1&&flag2==0)||(flag1==0&&flag2==1))
 	{
 		//find adjacent node to one not in graph which exsist in graph
	 	if(flag1==0)
	 	{
	 		k=0;
	 		int min=INT_MAX;
	 		int carnode;
	 		serachadjcentnode_car(s,graph);//carpathnodes[]
		 	dijkstra_dist(graph,s,'w',dist);
			 for(i=0;i<k;i++)
			 {
				 if(dist[pathnodes[i]]<min)
			 	{
			 		min=dist[pathnodes[i]];
			 		carnode=pathnodes[i];
			 	}
			 }

			 dijkstra_path(s,carnode,graph,dp,&dsize,'w');
			 int j=0;
			 for(i=0;i<dsize;i++)
				{
					path[j]=dp[i];
					j++;
				}
			dijkstra_path(carnode,d,graph,dp,&dsize,'c');
			for(i=1;i<dsize;i++)
			{
				path[j]=dp[i];
				j++;
			}
			*pathsize=j;
	 	}
	 	if(flag2==0)
	 	{	k=0;
	 		int min=INT_MAX;
	 		int carnode;
	 		serachadjcentnode_car(d,graph);
		 	dijkstra_dist(graph,d,'w',dist);
			 for(i=0;i<k;i++)
			 {
				 if(dist[pathnodes[i]]<min)
			 	{
			 		min=dist[pathnodes[i]];
			 		carnode=pathnodes[i];
			 	}
			 }

			 dijkstra_path(s,carnode,graph,dp,&dsize,'c');
			 int j=0;
			 for(i=0;i<dsize;i++)
				{
					path[j]=dp[i];
					j++;
				}
			dijkstra_path(carnode,d,graph,dp,&dsize,'w');
			for(i=1;i<dsize;i++)
			{
				path[j]=dp[i];
				j++;
			}
			*pathsize=j;
			cout<<"hi";
	 	}
 	}
}

int path_walk_public_walk(int s,int d,struct Graph* graph,int path[])
{
		k=0;
		struct qnode *walk_public,*public_public,*public_walk;
		struct qnode *min=intalize_min();
		int i;
	 	serachadjcentnode_public(s,graph);//publicpathnodes[]
		for(i=0;i<k;i++)
			 {
	 			walk_public=a_search(s,pathnodes[i],graph,'w');
			 	if(min->gprev>walk_public->gprev)
			 		{
			 			min=walk_public;
			 		}
			 }
			 walk_public=min;

		min=intalize_min();

		serachadjcentnode_public(d,graph);//publicpathnodes[]
		for(i=0;i<k;i++)
			 {
	 			public_walk=a_search(pathnodes[i],d,graph,'w');
			 	if(min->gprev>public_walk->gprev)
			 		{
			 			min=public_walk;
			 		}
			 }
			 public_walk=min;

		public_public=a_search(walk_public->path[walk_public->psize-1],public_walk->path[0],graph,'p');
		//printpath_asearch(walk_public->path,walk_public->psize,0);
		//cout<<"\nthe path by public:";
		//printpath_asearch(public_public->path,public_public->psize);
		//printpath_asearch(public_walk->path,public_walk->psize,0);
		int j;

		for( j=0;j<walk_public->psize;j++)
		{
			path[j]=walk_public->path[j];
		}

		for(i=1;i<public_public->psize;i++)
		{
			path[j]=public_public->path[i];
			j++;
		}

		for(i=1;i<public_walk->psize;i++)
		{
			path[j]=public_walk->path[i];
			j++;
		}
		return j;
}


void path_public_graph_a_search(int s,int d,struct Graph* graph,int path[],int *pathsize)
{
	struct qnode *p;
	struct qnode *min=(struct qnode*) malloc(sizeof(struct qnode));
	int flag1=checknodepresent_publicpath(s,graph);
	int flag2=checknodepresent_publicpath(d,graph);
 	int i;
 	min->gprev=range;
 	min->fvalue=range;
 	min->presentnode=range;
 	min->path[0]=range;
 	min->psize=range;

 	if(flag1==0&&flag2==0)
 	{
 		*pathsize=path_walk_public_walk(s,d,graph,path);
 	}
 	if(flag1==1&&flag2==1)
 	{
 		min=a_search(s,d,graph,'p');
 		for(i=0;i<min->psize;i++)
		{
			path[i]=min->path[i];
		}
		*pathsize=i;
 		//cout<<"\nthe path by public:";
 		//printpath_asearch(min->path,min->psize);
 	}

 	if((flag1==1&&flag2==0)||(flag1==0&&flag2==1))
 	{
 		//find adjacent node to one not in graph which exsist in graph
	 	if(flag1==0)
	 	{
	 		k=0;
	 		serachadjcentnode_public(s,graph);//publicpathnodes[]
	 		for(i=0;i<k;i++)
			 {
			 	p=a_search(s,pathnodes[i],graph,'w');

			 	if(min->gprev>p->gprev)
			 		{
			 			min=p;
			 		}
			 }
		//	printpath_asearch_public(min->path,min->psize);
			p=a_search(min->path[min->psize-1],d,graph,'p');
		//	cout<<"\nthe path by public";
		//	printpath_asearch(p->path,p->psize);
			int j=0;
			for(i=0;i<min->psize;i++)
			{
				path[j]=min->path[i];
				j++;
			}
				for(i=1;i<p->psize;i++)
			{
				path[j]=p->path[i];
				j++;
			}
			*pathsize=j;


	 	}
	 	if(flag2==0)
	 	{	k=0;
	 		serachadjcentnode_public(d,graph);
	 		for(i=0;i<k;i++)
			 {
	 			p=a_search(pathnodes[i],d,graph,'w');
			 	if(min->gprev>p->gprev)
			 		{
			 			min=p;
			 		}
			 }
			p=a_search(s,min->path[0],graph,'p');
//			cout<<"\nthe path by public:";
//			printpath_asearch(p->path,p->psize);
//			printpath_asearch(min->path,min->psize,0);
			int j=0;
			for(i=0;i<p->psize;i++)
			{
				path[j]=p->path[i];
				j++;
			}
			for(i=1;i<min->psize;i++)
				{
					path[j]=min->path[i];
					j++;
				}

			*pathsize=j;


	 	}

 	}
}

int path_walk_public_walk_d(int s,int d,struct Graph* graph,int path[])
{
		k=0;
		int i;
 		int min=INT_MAX;
		int dist[vertex];
 		int carnode,carnode2,dp[MAX],dsize,dp2[MAX],dsize2;
 		dijkstra_dist(graph,s,'w',dist);
		 serachadjcentnode_public(s,graph);//carpathnodes[]

		 for(i=0;i<k;i++)
		 {
			 if(dist[pathnodes[i]]<min)
		 	{
		 		min=dist[pathnodes[i]];
		 		carnode=pathnodes[i];
		 	}
		 }
	 	 dijkstra_path(s,carnode,graph,dp,&dsize,'w');

 		int j;
		for( j=0;j<dsize;j++)
		{
			path[j]=dp[j];
		}
		min=INT_MAX;
		serachadjcentnode_public(d,graph);//carpathnodes[]
	 	dijkstra_dist(graph,d,'w',dist);
		 for(i=0;i<k;i++)
		 {
			 if(dist[pathnodes[i]]<min)
		 	{
		 		min=dist[pathnodes[i]];
		 		carnode2=pathnodes[i];
		 	}
		 }
	 	 dijkstra_path(carnode2,d,graph,dp2,&dsize2,'w');
	 	 dijkstra_path(carnode,carnode2,graph,dp,&dsize,'p');
		cout<<"both done";
		for(i=1;i<dsize;i++)
		{
			path[j]=dp[i];
			j++;
		}

		for(i=1;i<dsize2;i++)
		{
			path[j]=dp2[i];
			j++;
		}
		return j;

}


void path_public_graph_d_search(int s,int d,struct Graph* graph,int path[],int *pathsize)
{
	int i;
	int flag1=checknodepresent_publicpath(s,graph);
	int flag2=checknodepresent_publicpath(d,graph);
	int dp[MAX],dsize;
	int dist[vertex];
 	if(flag1==0&&flag2==0)
 	{
        path_walk_public_walk_d(s,d,graph,path);
 	}
 	if(flag1==1&&flag2==1)
 	{
 		dijkstra_path(s,d,graph,path,&dsize,'p');
 		*pathsize=dsize;
 	}

 	if((flag1==1&&flag2==0)||(flag1==0&&flag2==1))
 	{
 		//find adjacent node to one not in graph which exsist in graph
	 	if(flag1==0)
	 	{
	 		k=0;
	 		int min=INT_MAX;
	 		int carnode;
	 		serachadjcentnode_public(s,graph);//carpathnodes[]
		 	dijkstra_dist(graph,s,'w',dist);
			 for(i=0;i<k;i++)
			 {
				 if(dist[pathnodes[i]]<min)
			 	{
			 		min=dist[pathnodes[i]];
			 		carnode=pathnodes[i];
			 	}
			 }

			 dijkstra_path(s,carnode,graph,dp,&dsize,'w');
			 int j=0;
			 for(i=0;i<dsize;i++)
				{
					path[j]=dp[i];
					j++;
				}
			dijkstra_path(carnode,d,graph,dp,&dsize,'p');
			for(i=1;i<dsize;i++)
			{
				path[j]=dp[i];
				j++;
			}
			*pathsize=j;
	 	}
	 	if(flag2==0)
	 	{	k=0;
	 		int min=INT_MAX;
	 		int carnode;
	 		serachadjcentnode_public(d,graph);
		 	dijkstra_dist(graph,d,'w',dist);
			 for(i=0;i<k;i++)
			 {
				 if(dist[pathnodes[i]]<min)
			 	{
			 		min=dist[pathnodes[i]];
			 		carnode=pathnodes[i];
			 	}
			 }

			 dijkstra_path(s,carnode,graph,dp,&dsize,'p');
			 int j=0;
			 for(i=0;i<dsize;i++)
				{
					path[j]=dp[i];
					j++;
				}
			dijkstra_path(carnode,d,graph,dp,&dsize,'w');
			for(i=1;i<dsize;i++)
			{
				path[j]=dp[i];
				j++;
			}
			*pathsize=j;
	 	}
 	}
}

void path_traffic_graph_d_search(int s,int d,struct Graph* graph,int path[],int *pathsize)
{
	int i;
	int flag1=checknodepresent_carpath(s,graph);
	int flag2=checknodepresent_carpath(d,graph);
	int dp[MAX],dsize;
	int dist[vertex];

 	if(flag1==1&&flag2==1)
 	{
 		dijkstra_path(s,d,graph,path,&dsize,'t');
 		*pathsize=dsize;
 	}
 	else
 	{
 		cout<<"\nNodes are not present in traffic path.";
 	}

}

void printpath_asearch(int arr[MAX],int size)
{

		for(int  i=0;i<size;i++)
		{
			cout<<"\t"<<arr[i];

 		}
}

void printpath_asearch(int arr[MAX],int size,int i)
{
 	cout<<"\nthe path by walk:";
		for(;i<size;i++)
		{
			cout<<"\t"<<arr[i];

 		}
}

void printpath_asearch_public(int arr[MAX],int size)
{
 	cout<<"\nthe path by walk:";
		for(int i=0;i<size;i++)
		{
			cout<<"\t"<<arr[i];

 		}
}

void print_final_path(int arr[],int size)
{
	cout<<"\nthe final path is:";
	for(int i=0;i<size;i++)
		{
			cout<<"\t"<<arr[i];

 		}
}

////////////////////////////////dijkstra/////////////////////////

//////////////////heap for dijkstra/////////////////////////////

struct dj_heap* newdj_heap(int v,int dist)
{
    struct dj_heap* min_heapNode =(struct dj_heap*) malloc(sizeof(struct dj_heap));
    min_heapNode->v = v;
    min_heapNode->dist = dist;
    return min_heapNode;
}

struct heap* createheap(int cap)
{
    struct heap* min_heap =(struct heap*) malloc(sizeof(struct heap));
    min_heap->pos = (int *)malloc(cap * sizeof(int));
    min_heap->size = 0;
    min_heap->cap = cap;
    min_heap->array =(struct dj_heap**) malloc(cap * sizeof(struct dj_heap*));

    return min_heap;
}

void swap_djnode(struct dj_heap** a, struct dj_heap** b)
{
    struct dj_heap* t = *a;
    *a = *b;
    *b = t;
}


void min_heapify(struct heap* min_heap, int i)
{
    int smallest, l, r;
    smallest = i;
    l  =left(i);
    r =right(i);

    if (l < min_heap->size &&
        min_heap->array[l]->dist < min_heap->array[smallest]->dist )
      smallest = l;

    if (r < min_heap->size &&
        min_heap->array[r]->dist < min_heap->array[smallest]->dist )
      smallest = r;

    if (smallest != i)
    {
        dj_heap *smallestNode = min_heap->array[smallest];
        dj_heap *iNode = min_heap->array[i];

        min_heap->pos[smallestNode->v] = i;
        min_heap->pos[iNode->v] = smallest;

        swap_djnode(&min_heap->array[smallest], &min_heap->array[i]);
        min_heapify(min_heap, smallest);
    }
}

int isEmpty(struct heap* min_heap)
{
    return min_heap->size == 1;
}
int isEmpty2(struct heap* min_heap)
{
    return min_heap->size == 2;
}

struct dj_heap* extract_Min(struct heap* min_heap)
{
    if (isEmpty(min_heap))
        return NULL;

    struct dj_heap* root = min_heap->array[0];

    struct dj_heap* lastNode = min_heap->array[min_heap->size - 1];
    min_heap->array[0] = lastNode;

    min_heap->pos[root->v] = min_heap->size-1;
    min_heap->pos[lastNode->v] = 0;

    --min_heap->size;
    min_heapify(min_heap, 0);

    return root;
}

void decrease_Key(struct heap* min_heap, int v, int dist)
{
    int i = min_heap->pos[v];

    min_heap->array[i]->dist = dist;

    while (i && min_heap->array[i]->dist < min_heap->array[parent(i)]->dist)
    {
        min_heap->pos[min_heap->array[i]->v] =parent(i);
        min_heap->pos[min_heap->array[parent(i)]->v] = i;
        swap_djnode(&min_heap->array[i],  &min_heap->array[parent(i)]);

        i = parent(i);
    }
}

bool isInheap(struct heap *min_heap, int v)
{
   if (min_heap->pos[v] < min_heap->size)
     return true;
   return false;
}

void dijkstra_dist(struct Graph* graph, int src,char type,int dist[vertex])
{
	//cout<<"in diji";
    int V = graph->V;
    struct heap* min_heap = createheap(V+1);
 	//V=V+1;
    for (int v=0;v<V;v++)
    {
        dist[v] = INT_MAX;
        min_heap->array[v]=newdj_heap(v,dist[v]);
        min_heap->pos[v]=v;
    }

    min_heap->array[src] = newdj_heap(src, dist[src]);
    min_heap->pos[src]   = src;
    dist[src] = 0;
    decrease_Key(min_heap, src, dist[src]);
 		    min_heap->size = V;


 	if(type=='c'||type=='C')
 	{int val=min_heap->size-cnode ;
  	 while (min_heap->size!= val)
  	// while (!isEmpty2(min_heap))
		    {

		        struct dj_heap* min_heapNode = extract_Min(min_heap);
		        int u = min_heapNode->v;
			   	serach_adjcent(u,graph,type);//change
			   	for(int i=0;i<sizeadjnodes;i++)
			   	  {
			   		int v=adjnodes[i];

			   		 int distance=dist[u]+car_weight(u,v,graph);//chaneg

		             if (isInheap(min_heap,v)&& dist[u] != INT_MAX&&distance < dist[v])
			            {
			                dist[v]=distance;
			                decrease_Key(min_heap,v,dist[v]);
			            }
			   	  }
		    }
 	}
 	if(type=='w'||type=='W')
 	{
			 while (!isEmpty(min_heap))
		    {

		        struct dj_heap* min_heapNode = extract_Min(min_heap);
		        int u = min_heapNode->v;
			   	serach_adjcent(u,graph,type);//change
			   	for(int i=0;i<sizeadjnodes;i++)
			   	  {
			   		int v=adjnodes[i];
				   int distance=dist[u]+walk_weight(u,v,graph);//chaneg

		             if (isInheap(min_heap,v)&& dist[u] != INT_MAX&&distance < dist[v])
			            {
			                dist[v]=distance;
			                decrease_Key(min_heap,v,dist[v]);
			            }
			   	  }
		    }
 	}
	if(type=='p'||type=='P')
 	{int val=min_heap->size-pnode ;
  	 while (min_heap->size!= val)
	  	   // while (!isEmpty(min_heap))
		    {
		        struct dj_heap* min_heapNode = extract_Min(min_heap);
		        int u = min_heapNode->v;
			   	serach_adjcent(u,graph,type);//change
			   	for(int i=0;i<sizeadjnodes;i++)
			   	  {
			   		int v=adjnodes[i];
			   		 int distance=dist[u]+public_weight(u,v,graph);//chaneg

		             if (isInheap(min_heap,v)&& dist[u] != INT_MAX&&distance < dist[v])
			            {
			                dist[v]=distance;
			                decrease_Key(min_heap,v,dist[v]);
			            }
			   	  }
		    }
 	}
	if(type=='d'||type=='D')
 	{// min_heap->size = V;
 		    while (!isEmpty(min_heap))
		    {
		        struct dj_heap* min_heapNode = extract_Min(min_heap);
		        int u = min_heapNode->v;
			   	serach_adjcent(u,graph,'w');//change
			   	for(int i=0;i<sizeadjnodes;i++)
			   	  {
			   		int v=adjnodes[i];
			   		 int distance=dist[u]+distance_weight(u,v,graph);//chaneg

		             if (isInheap(min_heap,v)&& dist[u] != INT_MAX&&distance < dist[v])
			            {
			                dist[v]=distance;
			                decrease_Key(min_heap,v,dist[v]);
			            }
			   	  }
		    }
 	}
  if(type=='T'||type=='t')
 	{int val=min_heap->size-cnode ;
  	 while (min_heap->size!= val)
  	 	// while (!isEmpty(min_heap))
		    {
		        struct dj_heap* min_heapNode = extract_Min(min_heap);
		        int u = min_heapNode->v;
			   	serach_adjcent(u,graph,'c');//change
			   	for(int i=0;i<sizeadjnodes;i++)
			   	  {
  		   		   		int v=adjnodes[i];
   					 	int distance=dist[u]+traffic_weight(u,v,graph);//change

				             if (isInheap(min_heap,v)&& dist[u] != INT_MAX&&distance < dist[v])
					            {
					                dist[v]=distance;
					                decrease_Key(min_heap,v,dist[v]);
		    			        }
			   	     }
			}
	return;
	}
}


void dijkstra(struct Graph* graph, int src,int prev[],char type)
{
	//cout<<"\nin diji";
    int V = graph->V;
   	//V=V+1;
    struct heap* min_heap = createheap(V);
 	int dist[vertex];
    for (int v=0;v<V;v++)
    {
        dist[v] = INT_MAX;
        min_heap->array[v]=newdj_heap(v,dist[v]);
        min_heap->pos[v]=v;

    }

    min_heap->array[src] = newdj_heap(src, dist[src]);
    min_heap->pos[src]   = src;
    dist[src] = 0;
    decrease_Key(min_heap, src, dist[src]);
	min_heap->size = V;


 	if(type=='c'||type=='C')
	{
//			cout<<"hi";
		    //	cout<<"cnode:  "<<cnode<<" c is: "<<c;
		    //	cout<<"min heap size"<<min_heap->size;
	int val=min_heap->size-cnode ;
 	while (min_heap->size!= val)
  	 //	 while (!isEmpty(min_heap))
		    {

		        struct dj_heap* min_heapNode = extract_Min(min_heap);
		        int u = min_heapNode->v;
			   	serach_adjcent(u,graph,type);//change
			   	for(int i=0;i<sizeadjnodes;i++)
			   	  {
  		   		   		int v=adjnodes[i];
   					 	int distance=dist[u]+car_weight(u,v,graph);//chaneg

				             if (isInheap(min_heap,v)&& dist[u] != INT_MAX&&distance < dist[v])
					            {
					                dist[v]=distance;
					                decrease_Key(min_heap,v,dist[v]);
					                prev[adjnodes[i]]=u;
		    			        }
			   	     }
			}//cout<<"  bye";
	return;
	}


 	if(type=='w'||type=='W')
 	{  //min_heap->size = V;
 		    while (!isEmpty(min_heap))
		    {
		        struct dj_heap* min_heapNode = extract_Min(min_heap);
		        int u = min_heapNode->v;
			   	serach_adjcent(u,graph,type);//change
			   	for(int i=0;i<sizeadjnodes;i++)
			   	  {

				   	int v=adjnodes[i];

			   		 int distance=dist[u]+walk_weight(u,v,graph);//chaneg

		             if (isInheap(min_heap,v)&& dist[u] != INT_MAX&&distance < dist[v])
			            {
			                dist[v]=distance;
			                decrease_Key(min_heap,v,dist[v]);
			                prev[adjnodes[i]]=u;
			            }
			   	  }
		    }
 	return;
 	}

	if(type=='p'||type=='P')
 	{
	 int val=min_heap->size-pnode ;
  	 while (min_heap->size!= val)
	  	    //while (!isEmpty(min_heap))
		    {
		        struct dj_heap* min_heapNode = extract_Min(min_heap);
		        int u = min_heapNode->v;
			   	serach_adjcent(u,graph,type);//change
			   	for(int i=0;i<sizeadjnodes;i++)
			   	  {
			   		int v=adjnodes[i];

			   		 int distance=dist[u]+public_weight(u,v,graph);//chaneg

		             if (isInheap(min_heap,v)&& dist[u] != INT_MAX&&distance < dist[v])
			            {
			                dist[v]=distance;
			                decrease_Key(min_heap,v,dist[v]);
			                prev[adjnodes[i]]=u;
			            }
			   	  }
		    }
		    //cout<<"bye";
		    return;
 	}
	if(type=='d'||type=='D')
 	{ //min_heap->size = V;
 		    while (!isEmpty(min_heap))
		    {

		        struct dj_heap* min_heapNode = extract_Min(min_heap);
		        int u = min_heapNode->v;
			   	serach_adjcent(u,graph,'w');//change
			   	for(int i=0;i<sizeadjnodes;i++)
			   	  {
			   		int v=adjnodes[i];

			   		 int distance=dist[u]+distance_weight(u,v,graph);//chaneg

		             if (isInheap(min_heap,v)&& dist[u] != INT_MAX&&distance < dist[v])
			            {
			                dist[v]=distance;
			                decrease_Key(min_heap,v,dist[v]);
			                prev[adjnodes[i]]=u;
			            }
			   	  }
		    }
		    return;
 	}
 	if(type=='T'||type=='t')
 	{
	 int val=min_heap->size-cnode ;
  	 while (min_heap->size!= val)
  	 	// while (!isEmpty(min_heap))
		    {
		        struct dj_heap* min_heapNode = extract_Min(min_heap);
		        int u = min_heapNode->v;
			   	serach_adjcent(u,graph,'c');//change
			   	for(int i=0;i<sizeadjnodes;i++)
			   	  {
  		   		   		int v=adjnodes[i];
   					 	int distance=dist[u]+traffic_weight(u,v,graph);//change

				             if (isInheap(min_heap,v)&& dist[u] != INT_MAX&&distance < dist[v])
					            {
					                dist[v]=distance;
					                decrease_Key(min_heap,v,dist[v]);
					                prev[adjnodes[i]]=u;
		    			        }
			   	     }
			}
	return;
	}
}

void  dijkstra_path(int src,int des,struct Graph *graph ,int dpath[],int *dsize,char type)
{

	 int path[MAX],prev[MAX];
	 dijkstra(graph,src,prev,type);
	 //cout<<"end diji";
	 int i=0;
	 int index=des;
	 while(index!=src)
	 {
	 	if(i==0)
	 	{
	 		path[i]=des;
	 		i++;
	 	}
	 	else
	 	{
	 		path[i]=prev[index];
	 		i++;
	 		index=prev[index];
	 	}
	 }
	 for(int j=0,k=i-1;j<i;j++,k--)
	 {
	 dpath[j]=path[k];
	 }
	*dsize=i;
}














//Initialize Vertex List i.e. assigns values to all vertex
void init_Vertex_List(Vertex_List *a)
{
    addVertex(a,1,"acacia",'L',233,114);
    addVertex(a,2,"abbey",'L',233,150);
    addVertex(a,3,"raine",'L',274,150);
    addVertex(a,4,"luxury",'H',233,164);
    addVertex(a,5,"ramac",'L',271,164);
    addVertex(a,6,"reighli",'L',272,192);
    addVertex(a,7,"kingston",'L',236,192);
    addVertex(a,8,"rabbit park",'L',236,230);
    addVertex(a,9,"keeble",'L',271,230);
    addVertex(a,10,"grand",'H',211,230);

    addVertex(a,11,"keel",'L',299,230);
    addVertex(a,12,"kent ",'L',298,183);
    addVertex(a,13,"kara",'L',346,183);
    addVertex(a,14,"drain",'L',346,227);
    addVertex(a,15,"oriental",'H',430,227);
    addVertex(a,16,"bounjour pp",'P',417,257);
    addVertex(a,17,"karen",'L',412,276);
    addVertex(a,18,"manchester",'L',394,274);
    addVertex(a,19,"merton",'L',272,273);
    addVertex(a,20,"harvich ",'P',235,273);

    addVertex(a,21,"brent",'L',211,273);
    addVertex(a,22,"hamlets",'L',254,290);
    addVertex(a,23,"hounslow",'L',254,324);
    addVertex(a,24,"west whites",'P',233,324);
    addVertex(a,25,"marriot",'H',273,324);
    addVertex(a,26,"camden",'L',298,354);
    addVertex(a,27,"ritz palace",'H',416,346);
    addVertex(a,28,"southampton",'L',460,355);
    addVertex(a,29,"wilton",'L',449,325);
    addVertex(a,30,"wardour",'L',462,301);

    addVertex(a,31,"warren",'L',435,314);
    addVertex(a,32,"whitefield",'L',414,313);
    addVertex(a,33,"victoria",'L',393,301);
    addVertex(a,34,"vere",'L',435,279);
    addVertex(a,35,"reserve",'H',470,273);
    addVertex(a,36,"paul farrow",'L',487,281);
    addVertex(a,37,"ian house",'L',471,253);
    addVertex(a,38,"byrne",'L',451,227);
    addVertex(a,39,"scott dean",'L',498,220);
    addVertex(a,40,"amoco",'L',534,277);

    addVertex(a,41,"tesco",'L',564,259);
    addVertex(a,42,"adams",'L',555,204);
    addVertex(a,43,"asda",'P',523,205);
    addVertex(a,44,"mortso",'L',412,374);
    addVertex(a,45,"saville",'L',432,377);
    addVertex(a,46,"regent",'L',168,162);
    addVertex(a,47,"queensway",'L',112,205);
    addVertex(a,48,"russel",'P',112,190);
    addVertex(a,49,"savoy palace",'L',100,190);
    addVertex(a,50,"smith square",'L',101,162);

    addVertex(a,51,"stevebiko",'L',119,143);
    addVertex(a,52,"villers",'L',112,178);
    addVertex(a,53,"vincent",'L',127,210);
    addVertex(a,54,"vinocia",'L',121,262);
    addVertex(a,55,"alton",'P',97,277);
    addVertex(a,56,"eelhi",'L',68,277);
    addVertex(a,57,"paerd",'L',36,282);
    addVertex(a,58,"portso",'L',64,219);
    addVertex(a,59,"ormond",'L',57,334);
    addVertex(a,60,"oxford",'L',62,362);

    addVertex(a,61,"villa",'H',107,358);
    addVertex(a,62,"pigot",'L',237,348);
    addVertex(a,63,"pleyden",'L',273,346);
    addVertex(a,64,"mathew",'L',274,366);
    addVertex(a,65,"malet",'L',273,401);
    addVertex(a,66,"grays",'P',294,392);
    addVertex(a,67,"west whales",'P',623,251);
    addVertex(a,68,"hamilton",'L',564,324);
    addVertex(a,69,"hardy",'L',452,110);
    addVertex(a,70,"apple green pp",'L',463,92);

    addVertex(a,71,"goodge",'L',439,73);
    addVertex(a,72,"florida",'L',378,73);
    addVertex(a,73,"garth",'L',331,65);
    addVertex(a,74,"frith",'L',321,88);
    addVertex(a,75,"drury",'L',332,87);
    addVertex(a,76,"darce",'L',340,101);
    addVertex(a,77,"bp connect ",'P',362,104);
    addVertex(a,78,"connaught",'L',365,125);
    addVertex(a,79,"caxton",'L',352,128);
    addVertex(a,80,"rowan",'L',352,141);

    addVertex(a,81,"northcout",'L',347,152);
    addVertex(a,82,"john street",'L',349,120);
    addVertex(a,83,"long lane",'L',339,109);
    addVertex(a,84,"baker",'L',339,119);
    addVertex(a,85,"bedfort",'L',325,110);
    addVertex(a,86,"bond",'L',112,276);
    addVertex(a,87,"cowcross",'L',96,244);
    addVertex(a,88,"halam",'L',101,228);
    addVertex(a,89,"bob marley",'L',84,222);
    addVertex(a,90," royal",'H',74,223);

    addVertex(a,91,"A4",'L',73,244);
    addVertex(a,92,"westin",'H',87,242);
    addVertex(a,93,"clerkenwell",'L',75,263);
    addVertex(a,94,"fann",'L',64,258);
    addVertex(a,95,"wihatterkar",'L',120,361);
    addVertex(a,96,"westcross",'L',132,370);
    addVertex(a,97,"camberwell",'L',119,368);
    addVertex(a,98,"wtling",'L',445,98);
    addVertex(a,99,"verginia",'L',436,99);
    addVertex(a,100,"tiara",'H',428,93);

    addVertex(a,101,"A21",'L',422,78);
    addVertex(a,102,"vegas",'L',492,230);
    addVertex(a,103,"clifford",'L',482,217);
    addVertex(a,104,"turnmill",'L',470,226);
    addVertex(a,105,"vigo",'L',471,240);
    addVertex(a,106,"golden",'L',481,226);
    addVertex(a,107,"hackney",'L',481,241);
    addVertex(a,108,"burlington",'L',449,192);
    addVertex(a,109,"boyle",'L',501,272);
    addVertex(a,110,"Beach",'L',429,244);

    addVertex(a,111,"A200",'L',447,245);

}

void Add_Roads(Graph* graph)
{
    addEdge(graph,96,97,"a",0.5*dist_factor,'W',1);
    addEdge(graph,97,95,"b",0.2*dist_factor,'W',1);
    addEdge(graph,95,61,"c",0.6*dist_factor,'W',1);
    addEdge(graph,61,60,"d",1.8*dist_factor,'P',1);
    addEdge(graph,61,86,"e",3.5*dist_factor,'P',1);
    addEdge(graph,60,59,"f",1*dist_factor,'P',1);
    addEdge(graph,59,57,"g",2.5*dist_factor,'P',1);
    addEdge(graph,57,56,"a",1.4*dist_factor,'P',1);
    addEdge(graph,56,55,"b",1*dist_factor,'P',1);
    addEdge(graph,55,86,"c",0.8*dist_factor,'P',1);
    addEdge(graph,86,54,"d",0.8*dist_factor,'P',1);
    addEdge(graph,57,58,"e",3.6*dist_factor,'W',1);
    addEdge(graph,56,93,"f",0.5*dist_factor,'W',1);
    addEdge(graph,55,87,"g",1.4*dist_factor,'P',1);
    addEdge(graph,93,94,"a",0.5*dist_factor,'W',1);
    addEdge(graph,93,92,"b",1.1*dist_factor,'W',1);
    addEdge(graph,94,91,"c",0.8*dist_factor,'W',1);
    addEdge(graph,91,92,"d",0.5*dist_factor,'W',1);
    addEdge(graph,92,87,"e",0.4*dist_factor,'W',1);
    addEdge(graph,91,90,"f",1*dist_factor,'W',1);
    addEdge(graph,92,89,"g",0.9*dist_factor,'W',1);
    addEdge(graph,87,88,"a",0.5*dist_factor,'W',1);
    addEdge(graph,90,89,"b",0.4*dist_factor,'W',1);
    addEdge(graph,89,88,"c",0.5*dist_factor,'W',1);
    addEdge(graph,54,53,"d",2.1*dist_factor,'P',1);
    addEdge(graph,53,47,"e",1*dist_factor,'P',2);
    addEdge(graph,88,47,"f",1*dist_factor,'W',1);
    addEdge(graph,47,48,"g",0.8*dist_factor,'P',2);
    addEdge(graph,48,49,"a",0.5*dist_factor,'P',2);
    addEdge(graph,48,52,"b",0.8*dist_factor,'C',1);
    addEdge(graph,49,50,"c",1.2*dist_factor,'P',2);
    addEdge(graph,50,51,"d",1.4*dist_factor,'C',1);
    addEdge(graph,50,46,"e",2.8*dist_factor,'P',2);

    addEdge(graph,65,66,"f",1*dist_factor,'C',3);
    addEdge(graph,65,64,"g",1.5*dist_factor,'C',3);
    addEdge(graph,66,26,"a",1.4*dist_factor,'C',2);
    addEdge(graph,26,63,"b",1.2*dist_factor,'P',2);
    addEdge(graph,64,63,"c",1*dist_factor,'C',3);
    addEdge(graph,62,24,"d",1*dist_factor,'C',1);
    addEdge(graph,63,25,"e",1.2*dist_factor,'P',1);
    addEdge(graph,25,23,"f",1*dist_factor,'P',2);
    addEdge(graph,24,23,"g",1*dist_factor,'C',2);
    addEdge(graph,23,22,"a",1.2*dist_factor,'P',3);
    addEdge(graph,22,20,"b",1*dist_factor,'P',3);
    addEdge(graph,20,19,"c",1.5*dist_factor,'C',1);
    addEdge(graph,20,21,"d",1*dist_factor,'C',2);
    addEdge(graph,21,10,"e",1.7*dist_factor,'C',1);
    addEdge(graph,20,8,"f",1.7*dist_factor,'P',1);
    addEdge(graph,10,8,"g",1*dist_factor,'C',1);
    addEdge(graph,8,9,"a",1.5*dist_factor,'P',1);
    addEdge(graph,9,11,"b",1.2*dist_factor,'P',2);
    addEdge(graph,11,12,"c",1.7*dist_factor,'P',2);
    addEdge(graph,9,6,"d",1.5*dist_factor,'P',3);
    addEdge(graph,8,7,"e",1.5*dist_factor,'C',1);
    addEdge(graph,7,6,"f",1.5*dist_factor,'C',2);
    addEdge(graph,6,5,"g",1.2*dist_factor,'P',3);
    addEdge(graph,5,3,"a",0.3*dist_factor,'C',3);
    addEdge(graph,5,4,"b",1.5*dist_factor,'P',2);
    addEdge(graph,4,2,"c",0.7*dist_factor,'P',2);
    addEdge(graph,2,3,"d",1.6*dist_factor,'C',2);
    addEdge(graph,2,1,"e",1.3*dist_factor,'P',2);

    addEdge(graph,44,45,"f",1*dist_factor,'W',1);
    addEdge(graph,44,27,"g",1.1*dist_factor,'W',1);
    addEdge(graph,45,28,"a",1.3*dist_factor,'W',1);
    addEdge(graph,27,28,"b",2*dist_factor,'W',1);
    addEdge(graph,28,29,"c",1.5*dist_factor,'W',1);
    addEdge(graph,27,32,"d",1.4*dist_factor,'C',1);
    addEdge(graph,32,31,"e",1*dist_factor,'C',2);
    addEdge(graph,32,33,"f",1*dist_factor,'C',2);
    addEdge(graph,29,30,"g",1*dist_factor,'P',2);
    addEdge(graph,30,36,"a",1.2*dist_factor,'P',1);
    addEdge(graph,33,18,"b",1*dist_factor,'C',1);
    addEdge(graph,32,17,"c",1.2*dist_factor,'C',1);
    addEdge(graph,31,34,"d",1.2*dist_factor,'C',3);
    addEdge(graph,18,17,"e",1*dist_factor,'C',2);
    addEdge(graph,34,35,"f",1.3*dist_factor,'C',3);
    addEdge(graph,35,36,"g",1*dist_factor,'C',3);
    addEdge(graph,36,40,"a",2.2*dist_factor,'P',1);
    addEdge(graph,40,109,"b",1.2*dist_factor,'P',1);
    addEdge(graph,17,16,"c",1*dist_factor,'C',3);
    addEdge(graph,40,41,"d",1.2*dist_factor,'P',1);
    addEdge(graph,35,37,"e",0.8*dist_factor,'C',2);
    addEdge(graph,16,15,"f",1.2*dist_factor,'C',1);
    addEdge(graph,37,38,"g",1.2*dist_factor,'C',1);
    addEdge(graph,15,110,"a",0.8*dist_factor,'C',1);
    addEdge(graph,110,111,"b",0.8*dist_factor,'W',1);
    addEdge(graph,37,105,"c",0.5*dist_factor,'C',1);
    addEdge(graph,105,107,"d",0.4*dist_factor,'W',1);
    addEdge(graph,105,104,"e",0.4*dist_factor,'W',1);
    addEdge(graph,107,106,"f",0.4*dist_factor,'W',1);
    addEdge(graph,106,102,"g",0.4*dist_factor,'W',1);
    addEdge(graph,106,104,"a",0.4*dist_factor,'W',1);
    addEdge(graph,106,103,"b",0.2*dist_factor,'W',1);
    addEdge(graph,103,39,"c",0.6*dist_factor,'W',1);
    addEdge(graph,109,39,"d",2.2*dist_factor,'C',1);
    addEdge(graph,40,43,"e",2.6*dist_factor,'P',1);
    addEdge(graph,41,43,"f",2.8*dist_factor,'P',1);
    addEdge(graph,41,42,"g",2.4*dist_factor,'P',1);
    addEdge(graph,15,108,"a",1.4*dist_factor,'P',1);
    addEdge(graph,108,43,"b",3.6*dist_factor,'P',1);
    addEdge(graph,42,67,"c",3.4*dist_factor,'P',1);
    addEdge(graph,67,68,"d",3.5*dist_factor,'P',1);

    addEdge(graph,81,80,"e",0.5*dist_factor,'W',1);
    addEdge(graph,80,79,"f",0.5*dist_factor,'W',1);
    addEdge(graph,79,78,"g",0.5*dist_factor,'W',1);
    addEdge(graph,79,82,"a",0.3*dist_factor,'W',1);
    addEdge(graph,82,84,"b",0.2*dist_factor,'W',1);
    addEdge(graph,84,83,"c",0.3*dist_factor,'W',1);
    addEdge(graph,83,85,"d",0.6*dist_factor,'W',1);
    addEdge(graph,83,76,"e",0.2*dist_factor,'W',1);
    addEdge(graph,78,77,"f",1.1*dist_factor,'W',1);
    addEdge(graph,77,76,"g",1*dist_factor,'W',1);
    addEdge(graph,76,75,"a",0.8*dist_factor,'C',1);
    addEdge(graph,75,74,"b",0.5*dist_factor,'P',1);
    addEdge(graph,75,73,"c",1.1*dist_factor,'P',2);
    addEdge(graph,73,72,"d",1.9*dist_factor,'P',2);
    addEdge(graph,72,71,"e",2.4*dist_factor,'P',2);
    addEdge(graph,71,70,"f",1.2*dist_factor,'P',2);
    addEdge(graph,70,69,"g",1*dist_factor,'P',2);
    addEdge(graph,69,98,"a",0.8*dist_factor,'C',1);
    addEdge(graph,98,99,"b",0.3*dist_factor,'W',1);
    addEdge(graph,99,100,"c",0.2*dist_factor,'W',1);
    addEdge(graph,100,101,"d",1*dist_factor,'W',1);

    addEdge(graph,69,108,"e",3.4*dist_factor,'P',1);
    addEdge(graph,1,74,"f",3.9*dist_factor,'P',1);
    addEdge(graph,12,13,"g",2*dist_factor,'P',1);
    addEdge(graph,13,14,"a",1.4*dist_factor,'P',1);
    addEdge(graph,14,15,"b",3.5*dist_factor,'P',1);
    addEdge(graph,27,26,"c",5*dist_factor,'C',1);
    addEdge(graph,46,4,"d",2.8*dist_factor,'P',1);

}



