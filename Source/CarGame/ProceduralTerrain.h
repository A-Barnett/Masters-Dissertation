#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "ProceduralMeshComponent.h"
#include "Async/Async.h"
#include "Containers/Queue.h"
#include "Misc/ScopeLock.h"
#include "Templates/Function.h" 
#include "HAL/CriticalSection.h"
#include "RealtimeMeshActor.h"
#include "GameFramework/Actor.h"
#include "RealtimeMeshComponent.h"
#include "RealtimeMeshSimple.h"
#include "Misc/ScopeLock.h"
#include "Templates/Function.h"
#include "CoreMinimal.h"
#include "RealtimeMeshCore.h"
#include "RealtimeMeshConfig.h"




#include "ProceduralTerrain.generated.h"

class TerrainComponent {
private:
    UProceduralMeshComponent* MeshComponent;
    FVector2D GridPosition;
    int LOD;
    bool IsActive = true;
    int Index;
    bool IsInitialised = false;
    FRealtimeMeshSectionGroupKey GroupKey;
    FRealtimeMeshSectionKey Key;

public:
    TerrainComponent(UProceduralMeshComponent* InMeshComponent, FVector2D InGridPosition, int InLOD, int Inindex)
        : MeshComponent(InMeshComponent), GridPosition(InGridPosition), LOD(InLOD), IsActive(true), Index(Inindex),IsInitialised(false) {
    }
    UProceduralMeshComponent* GetMeshComponent() const {
        return MeshComponent;
    }
    void SetMeshComponent(UProceduralMeshComponent* InMeshComponent) {
        MeshComponent = InMeshComponent;
    }
    FVector2D GetGridPosition() const {
        return GridPosition;
    }
    void SetGridPosition(const FVector2D& InGridPosition) {
        GridPosition = InGridPosition;
    }
    int GetLOD() const {
        return LOD;
    }
    void SetLOD(int InLOD) {
        LOD = InLOD;
    }
    bool GetIsActive() const {
        return IsActive;
    }
    void SetIsActive(bool InIsActive) {
        IsActive = InIsActive;
    }
    int GetIndex() const {
        return Index;
    }
    void SetIndex(int InIndex) {
        Index = InIndex;
    }
    bool GetIsInitialised() const {
        return IsInitialised;
    }
    void SetIsInitialised(bool InIsInitialised) {
        IsInitialised = InIsInitialised;
    }
    const FRealtimeMeshSectionGroupKey GetGroupKey() {
        return GroupKey;
    }
    void SetGroupKey(const FRealtimeMeshSectionGroupKey InGroupKey) {
        GroupKey = InGroupKey;
    }

    const FRealtimeMeshSectionKey GetKey() {
        return Key;
    }
    void SetKey(const FRealtimeMeshSectionKey InKey) {
        Key = InKey;
    }
};

class MeshGenerationFactory
{
private:
    TQueue<TFunction<void()>> MeshGenerationQueue; // Queue of mesh creation tasks
    FCriticalSection Mutex;                        // Mutex to ensure thread safety

public:
    void EnqueueMeshTask(TFunction<void()> MeshTask)
    {
        FScopeLock Lock(&Mutex);
        MeshGenerationQueue.Enqueue(MoveTemp(MeshTask));
    }

    void ProcessMeshTasks()
    {
        TFunction<void()> MeshTask;
        if (MeshGenerationQueue.Dequeue(MeshTask))
        {
            MeshTask(); // Execute one task per frame
        }
    }

    bool HasPendingTasks() const
    {
        return !MeshGenerationQueue.IsEmpty();
    }
};




USTRUCT(BlueprintType)
struct FTerrainLayer
{
    GENERATED_BODY()

    UPROPERTY(EditAnywhere, Category = "Layer")
    float Amplitude;

    UPROPERTY(EditAnywhere, Category = "Layer")
    float Frequency;

    UPROPERTY(EditAnywhere, Category = "Layer")
    float Contribution;

    UPROPERTY(EditAnywhere, Category = "Layer")
    bool Enabled = true;
};


UCLASS()
class CARGAME_API AProceduralTerrain : public ARealtimeMeshActor
{
    GENERATED_BODY()

public:
    AProceduralTerrain();
    
protected:
    virtual void BeginPlay() override;

public:
    virtual void Tick(float DeltaTime) override;
   
    void SmoothPathPointsHeight(float smoothLevel);
    void UpdateTerrain();
    TerrainComponent* FindTerrainComponent(const FVector2D& GridPosition);
    TerrainComponent* CreateTerrainComponent(const FVector2D& GridPosition, int LOD);
    void RemoveTerrainComponent(TerrainComponent* terrain);
    void ResetPlayerPhysics();

    URealtimeMeshSimple* RealtimeMesh;
    int totalCreated = 0;


    MeshGenerationFactory TerrainMeshFactory;
    TPair<float, float> randomOffset;
    FVector2D PlayerGridPos = FVector2D(0,0);
    APawn* PlayerPawn;
    float currentGridSize;

    const float GridSize = 50.0f;
    TMap<FIntPoint, TArray<FVector>> PathGrid = {
    { FIntPoint(0, 0), { FVector(0, 0, 0) } }
    };

    UPROPERTY(EditAnywhere, Instanced, Category = "Terrain")
    TArray<UProceduralMeshComponent*> MeshComponent;
    TArray<TerrainComponent*> TerrainComponents;
    
    FCriticalSection Mutex;
    FCriticalSection ComponentMutex;
    FCriticalSection DataMutex;

    UPROPERTY(EditAnywhere, Category = "Terrain")
    int32 Width;

    UPROPERTY(EditAnywhere, Category = "Terrain")
    int32 Height = 100;

    UPROPERTY(EditAnywhere, Category = "Terrain")
    float Scale;

    UPROPERTY(EditAnywhere, Category = "Terrain")
    float UVScale = 1.0f;

    UPROPERTY(EditAnywhere, Category = "Terrain Layers")
    TArray<FTerrainLayer> TerrainLayers;

    UPROPERTY(EditAnywhere, Category = "Materials")
    UMaterialInterface* TerrainMaterial;

    UPROPERTY(EditAnywhere, Category = "Materials")
    UMaterialInterface* PathMaterial;

    UPROPERTY(EditAnywhere, Category = "Path")
    int32 NumPoints;

    UPROPERTY(EditAnywhere, Category = "Path")
    float Variance;

    UPROPERTY(EditAnywhere, Category = "Path")
    float Thickness;

    UPROPERTY(EditAnywhere, Category = "Path")
    int32 PathSeed; 

    UPROPERTY(EditAnywhere, Category = "Path")
    float HeightAdjust;

    UPROPERTY(EditAnywhere, Category = "Path")
    int32 ThicknessDetail = 2;

    UPROPERTY(EditAnywhere, Category = "Path")
    float PathTextureScale = 1;

    UPROPERTY(EditAnywhere, Category = "Path")
    float EdgeHeightOffset = -1;

    UPROPERTY(EditAnywhere, Category = "Path")
    float Flatness = 0;

    UPROPERTY(EditAnywhere, Category = "Path")
    float ThicknessOffset = 0;

    UPROPERTY(EditAnywhere, Category = "Path")
    float SmoothingThicknessOffset = 0;

    UPROPERTY(EditAnywhere, Category = "Path")
    float PathHeightSmooth = 0;

    UPROPERTY(EditAnywhere, Category = "Path")
    float SmoothingSize = 12;

    UPROPERTY(EditAnywhere, Category = "Path")
    int TurnSmoothPoints = 20; // How many points to spread a turn across

    UPROPERTY(EditAnywhere, Category = "Path")
    int TurnSize = 100;// How often points to start turn

    UPROPERTY(EditAnywhere, Category = "Path")
    float SmoothingPasses = 0.5f;
    UFUNCTION(CallInEditor, Category = "Path")
    void RemoveTerrain()
    {
        randomOffset.Value = rand() % 100000;
        randomOffset.Key = rand() % 100000;

        // Array to store components to remove
        TArray<UProceduralMeshComponent*> ProceduralMeshComponents;

        // Iterate through all components of the actor
        TArray<UActorComponent*> Components;
        GetComponents(Components);

        for (UActorComponent* Component : Components)
        {
            if (UProceduralMeshComponent* ProcMesh = Cast<UProceduralMeshComponent>(Component))
            {
                ProceduralMeshComponents.Add(ProcMesh);
            }
        }

        // Remove all procedural mesh components
        for (UProceduralMeshComponent* ProcMesh : ProceduralMeshComponents)
        {
            ProcMesh->DestroyComponent();
        }

        // Clear stored mesh references
        MeshComponent.Empty();

        UE_LOG(LogTemp, Log, TEXT("All terrain procedural meshes removed."));
    }

    TArray<FVector> PathPoints;
    TArray<FVector> PathVertices;
    TArray<int32> PathTriangles;
    TArray<FVector> PathNormals;
    TArray<FVector2D> PathUVs;
    void GenerateTerrain();

private:
    float CalculateHeight(int32 X, int32 Y) const;
    float CalculateHeightOnPath(int32 X, int32 Y) const;
    float CalculateNoiseAtPoint(int32 X, int32 Y) const;
    void GeneratePath();
    bool IsOnPath(int32 X, int32 Y, bool useOffset) const;
    float DistFromPath(int32 X, int32 Y, bool useOffset) const;
    void GenerateTerrainSection(TerrainComponent* component);
    void GeneratePathMesh();
    void DisplayPathMesh();
};
