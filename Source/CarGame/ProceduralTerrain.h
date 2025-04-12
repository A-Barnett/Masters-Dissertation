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
#include <mutex>


#include "ProceduralTerrain.generated.h"

class TerrainComponent {
private:
    FVector2D GridPosition;
    int LOD;
    int OldLOD;
    bool IsActive = true;
    int32 Index;
    bool IsInitialised = false;
    FRealtimeMeshSectionGroupKey GroupKey;
    FRealtimeMeshSectionKey Key;
    TSharedPtr<TRealtimeMeshBuilderLocal<uint16, FPackedNormal, FVector2DHalf, 1>> PrevBuilder;
    TSharedPtr<FRealtimeMeshStreamSet> PrevStreamSet;
    bool isNewPos = true;
    std::mutex MeshMutex;

public:
    TerrainComponent( FVector2D InGridPosition, int InLOD, int32 Inindex)
        :  GridPosition(InGridPosition), LOD(InLOD), IsActive(true), Index(Inindex),IsInitialised(false) {
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
    int32 GetIndex() const {
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

    void SetPrevBuilder(TSharedPtr<TRealtimeMeshBuilderLocal<uint16, FPackedNormal, FVector2DHalf, 1>>&& InBuilder) {
        PrevBuilder = MoveTemp(InBuilder);
    }

    TSharedPtr<TRealtimeMeshBuilderLocal<uint16, FPackedNormal, FVector2DHalf, 1>> GetPrevBuilder() {
        return PrevBuilder;
    }

    void SetPrevStreamset(TSharedPtr<FRealtimeMeshStreamSet> StreamSetIn) {
        PrevStreamSet = MoveTemp(StreamSetIn);
    }

    void SetIsNewPos(bool isNewPosIn) {
        isNewPos = isNewPosIn;
    }

    bool GetIsNewPos() {
        return isNewPos;
    }

    int GetOldLOD() const {
        return OldLOD;
    }
    void SetOldLOD(int InLOD) {
        OldLOD = InLOD;
    }
    bool TryLockMesh() {
        return MeshMutex.try_lock();
    }
    void UnlockMesh() {
        MeshMutex.unlock();
    }

};

class PathComponent {
private:
    FRealtimeMeshSectionGroupKey GroupKey;
    FRealtimeMeshSectionKey Key;
    int32 Index;

public:
    PathComponent(int32 IndexIn) {
        Index = IndexIn;
    };
  
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
    int32 GetIndex() const {
        return Index;
    }
    void SetIndex(int InIndex) {
        Index = InIndex;
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
   
    void SmoothPathPointsHeight(float smoothLevel, int32 offset);
    void UpdateTerrain(bool initial = false);
    TerrainComponent* FindTerrainComponent(const FVector2D& GridPosition);
    TerrainComponent* CreateTerrainComponent(const FVector2D& GridPosition, int LOD);
    void ResetPlayerPhysics();
    void PlayerStartPos();

    UPROPERTY()
    URealtimeMeshComponent* TerrainMeshComponent;

    UPROPERTY()
    URealtimeMeshComponent* PathMeshComponent;

    URealtimeMeshSimple* RealtimeMesh;
    URealtimeMeshSimple* PathRealtimeMesh;

    int totalCreated = 0;

    TPair<float, float> randomOffset;
    FVector2D PlayerGridPos = FVector2D(0,0);
    APawn* PlayerPawn;
    float currentGridSize;

    const float GridSize = 50.0f;
    TMap<FIntPoint, TArray<FVector>> PathGrid = {
    { FIntPoint(0, 0), { FVector(0, 0, 0) } }
    };


    TArray<TerrainComponent*> TerrainComponents;
    TArray<PathComponent*> PathComponents;
    
    FCriticalSection Mutex;
    FCriticalSection ComponentMutex;
    FCriticalSection DataMutex;

    UPROPERTY(EditAnywhere, Category = "Terrain")
    int32 Width;

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

   

    TArray<FVector3f> PathPoints;
    TMap<FIntPoint, TArray<FVector3f>> PathGridMap;
    TMap<FIntPoint, TArray<FVector3f>> PathVertMap;
    void GenerateTerrain();

private:
    FVector2D CurrentPosition;
    FVector2D CurrentDirection;
    float CurrentTurnAngle;
    FRandomStream RandomStream;
    TArray<FVector3f> LastVertexes;
    int32 pointCount = 0;
    int32 pointsGenerated = 0;
    int32 makePerUpdate = 100;
    int32 prevEnd = 0;
    int32 vertsInPoint = 0;
    int32 pathUpdateCount = 0;
    int32 pathComponentsCount = 0;
    float CalculateHeight(int32 X, int32 Y, TArray<FVector3f>* PointsMap, TArray<FVector3f>* VertMap) const;
    float CalculateHeightOnPath(int32 X, int32 Y, TArray<FVector3f>* VertMap) const;
    float CalculateNoiseAtPoint(int32 X, int32 Y) const;
    void GeneratePath(bool start, int32 offset);
    bool IsOnPath(int32 X, int32 Y, bool useOffset, TArray<FVector3f>* PointsMap) const;
    float DistFromPath(int32 X, int32 Y, bool useOffset, TArray<FVector3f>* PointsMap) const;
    void GenerateTerrainSection(TerrainComponent* component);
    void GeneratePathMesh(int32 offset, PathComponent* path);
    void RemovePathMesh(PathComponent* path);
};
